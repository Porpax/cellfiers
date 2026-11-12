#' Plot neighbour-group comparisons with explicit test labels
#'
#' Plots the data used by [compare_neighbors()] and annotates each cell
#' type with its P value. The subtitle records the analysis and permutation
#' units so exploratory FOV-level results cannot be mistaken for clustered or
#' sample-level inference.
#'
#' @param data Grouped long FOV-level neighbour summary.
#' @param statistics Output of [compare_neighbors()].
#' @param value_col,group_col,sample_col,celltype_col Column names in `data`.
#' @param p_column P-value column in `statistics`. Defaults to the raw
#'   `"p.value"`; use `"p.adjust"` to display multiplicity-adjusted values.
#' @param percent Multiply proportions by 100 for display.
#' @param palette Named group colours.
#' @param geom One of `"both"`, `"violin"`, or `"boxplot"`.
#' @param show_points Show individual sample/FOV observations.
#' @param point_size,alpha Point size and transparency. Points use a black
#'   outline and the group colour as fill so they remain visible over violins.
#' @param violin_alpha Transparency of the violin fill.
#' @param box_color Colour of boxplot outlines, whiskers and median lines.
#' @param visual_iqr_filter If `TRUE`, observations outside the pooled
#'   cell-type-specific `Q1 - visual_iqr_multiplier * IQR` to
#'   `Q3 + visual_iqr_multiplier * IQR` interval are omitted from the figure
#'   only. Statistics supplied in `statistics` are not recalculated.
#' @param visual_iqr_multiplier IQR multiplier used only when
#'   `visual_iqr_filter = TRUE`.
#' @param title Optional title.
#' @param y_lab Optional y-axis label.
#'
#' @return A [ggplot2::ggplot()] object.
#' @examples
#' d <- expand.grid(Sample = paste0("S", 1:6), FOV = c("F1", "F2"),
#'                  CellType = "CD8 T", stringsAsFactors = FALSE)
#' d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
#' d$mean_prop <- seq_len(nrow(d)) / nrow(d)
#' tests <- compare_neighbors(d, analysis_level = "sample")
#' plot_neighbors(d, tests)
#' @export
plot_neighbors <- function(
    data,
    statistics,
    value_col = "mean_prop",
    group_col = "Group",
    sample_col = "Sample",
    celltype_col = "CellType",
    p_column = "p.value",
    percent = TRUE,
    palette = c(R = "#E64B35", NR = "#4DBBD5"),
    geom = c("both", "violin", "boxplot"),
    show_points = TRUE,
    point_size = 1.4,
    alpha = 0.65,
    violin_alpha = 0.4,
    box_color = "black",
    visual_iqr_filter = FALSE,
    visual_iqr_multiplier = 1.5,
    title = NULL,
    y_lab = NULL) {
  geom <- match.arg(geom)
  vc <- .find_col(data, value_col)
  gc <- .find_col(data, group_col)
  sc <- .find_col(data, sample_col)
  cc <- .find_col(data, celltype_col)
  if (any(vapply(list(vc, gc, sc, cc), is.null, logical(1)))) {
    stop("Could not find value, group, sample, and/or cell-type columns.",
         call. = FALSE)
  }
  required_stats <- c(
    "CellType", "analysis_level", "method", "method_label", p_column
  )
  if (!is.data.frame(statistics) || !all(required_stats %in% names(statistics))) {
    stop("`statistics` must be output from compare_neighbors().",
         call. = FALSE)
  }
  methods <- unique(statistics$method)
  levels <- unique(statistics$analysis_level)
  if (length(methods) != 1L || length(levels) != 1L) {
    stop("`statistics` must contain exactly one analysis method and level.",
         call. = FALSE)
  }

  d <- data.frame(
    Sample = as.character(data[[sc]]),
    Group = as.character(data[[gc]]),
    CellType = as.character(data[[cc]]),
    value = data[[vc]], stringsAsFactors = FALSE
  )
  d <- d[is.finite(d$value), , drop = FALSE]
  if (methods == "sample_wilcox") {
    d <- d |>
      dplyr::group_by(.data$Sample, .data$Group, .data$CellType) |>
      dplyr::summarise(value = mean(.data$value), .groups = "drop")
  }
  multiplier <- if (percent) 100 else 1
  d$value_plot <- d$value * multiplier
  group_order <- unique(c(statistics$group1, statistics$group2))
  d$Group <- factor(d$Group, levels = group_order)

  plot_d <- d
  if (isTRUE(visual_iqr_filter)) {
    if (!is.numeric(visual_iqr_multiplier) ||
        length(visual_iqr_multiplier) != 1L ||
        !is.finite(visual_iqr_multiplier) || visual_iqr_multiplier <= 0) {
      stop("`visual_iqr_multiplier` must be one positive finite number.",
           call. = FALSE)
    }
    plot_d <- plot_d |>
      dplyr::group_by(.data$CellType) |>
      dplyr::mutate(
        .__q1 = stats::quantile(.data$value_plot, 0.25, na.rm = TRUE),
        .__q3 = stats::quantile(.data$value_plot, 0.75, na.rm = TRUE),
        .__iqr = .data$.__q3 - .data$.__q1,
        .__keep = ifelse(
          is.finite(.data$.__iqr) & .data$.__iqr > 0,
          .data$value_plot >= .data$.__q1 - visual_iqr_multiplier * .data$.__iqr &
            .data$value_plot <= .data$.__q3 + visual_iqr_multiplier * .data$.__iqr,
          TRUE
        )
      ) |>
      dplyr::filter(.data$.__keep) |>
      dplyr::ungroup() |>
      dplyr::select(-dplyr::all_of(c(".__q1", ".__q3", ".__iqr", ".__keep")))
  }

  annotation <- data.frame(
    CellType = statistics$CellType,
    p = statistics[[p_column]], stringsAsFactors = FALSE
  )
  maxima <- plot_d |>
    dplyr::group_by(.data$CellType) |>
    dplyr::summarise(y = max(.data$value_plot, na.rm = TRUE), .groups = "drop")
  annotation <- dplyr::left_join(annotation, maxima, by = "CellType")
  annotation$y <- ifelse(annotation$y > 0, annotation$y * 1.10, 0.05)
  prefix <- if (p_column == "p.adjust") "BH P = " else "P = "
  annotation$label <- paste0(
    prefix, format.pval(annotation$p, digits = 2, eps = 2.2e-16)
  )

  method_subtitle <- unique(statistics$method_label)
  if (p_column == "p.adjust") {
    method_subtitle <- paste0(
      method_subtitle, "; multiple testing: ",
      unique(statistics$p_adjust_method)
    )
  }
  if (isTRUE(visual_iqr_filter)) {
    method_subtitle <- paste0(
      method_subtitle,
      "; display zoom: pooled cell-type IQR range (statistics use all data)"
    )
  }
  if (is.null(y_lab)) {
    y_lab <- if (percent) "% of centre-cell neighbours" else
      "Neighbour proportion"
  }

  p <- ggplot2::ggplot(
    plot_d,
    ggplot2::aes(
      x = .data$Group, y = .data$value_plot,
      color = .data$Group, fill = .data$Group
    )
  )
  if (geom %in% c("both", "violin")) {
    p <- p + ggplot2::geom_violin(alpha = violin_alpha, trim = TRUE)
  }
  if (geom %in% c("both", "boxplot")) {
    p <- p + ggplot2::geom_boxplot(
      width = 0.2, color = box_color, alpha = 0.35,
      linewidth = 0.6, outlier.shape = NA
    )
  }
  if (show_points) {
    p <- p + ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.08, height = 0),
      shape = 21, color = "black", stroke = 0.35,
      size = point_size, alpha = alpha
    )
  }
  p +
    ggplot2::geom_text(
      data = annotation,
      ggplot2::aes(x = 1.5, y = .data$y, label = .data$label),
      inherit.aes = FALSE, size = 3.1, color = "black"
    ) +
    ggplot2::facet_wrap(~CellType, nrow = 1, scales = "free_y") +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
    ggthemes::theme_few() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(fill = NA, color = "black"),
      legend.position = "none",
      strip.text = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(
      x = NULL, y = y_lab, title = title, subtitle = method_subtitle
    )
}

#' @rdname plot_neighbors
#' @export
plot_neighbor_comparison <- plot_neighbors
