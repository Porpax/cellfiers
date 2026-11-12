# encoding: UTF-8
#' Group comparison violin/box plot (original workflow style)
#'
#' Reproduces the group-comparison figure of the original workflow
#' (`mIF_Tzone_05_QC_Visual.R`): opaque violin plots filled with the group
#' colours, white boxplots (`width = 0.2, alpha = 0.5`), `ggthemes::theme_few`
#' with a black panel border, one facet row, a percentage axis
#' (`values x 100`, breaks every 2.5) and pairwise Wilcoxon p-values drawn by
#' `ggpubr::stat_compare_means` (`label = "p.format"`, `size = 3.5`).
#'
#' @param data A data frame with a numeric value column, a group column and
#'   optionally a `by` column.
#' @param value_col Name of the numeric column to plot.
#' @param group_col Name of the group column (x axis).
#' @param by Optional column whose levels are faceted in one row (e.g.
#'   `"CellType"`).
#' @param method Test used for p-values (passed to
#'   `ggpubr::stat_compare_means`), e.g. `"wilcox.test"` or `"t.test"`.
#' @param comparisons Optional list of group pairs to test, e.g.
#'   `list(c("R", "NR"))`. If `NULL`, all pairwise combinations are used.
#' @param geom Plot geometry: `"both"` (violin + box, default, as in the
#'   original), `"box"` or `"violin"`.
#' @param percent Logical; if `TRUE` (default, original behaviour), values are
#'   multiplied by 100 and the y axis uses breaks every 2.5.
#' @param y_lab Y axis label (default `"% of DC 20μm neighbors"` as in the
#'   original).
#' @param palette Optional character vector of colours, one per group. Use a
#'   named vector to fix colours to groups, e.g.
#'   `c(R = "#E64B35", NR = "#4DBBD5")`. Defaults to the NPG red/blue pair
#'   for two groups.
#' @param jitter Logical; if `TRUE`, individual points are shown (the original
#'   figure does not show them).
#' @param point_size,alpha Jitter point size and transparency.
#' @param title Optional plot title. The original figure has no title.
#' @param iqr_filter Logical; apply the original per-cell-type 1.5-IQR filter
#'   before plotting.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#'
#' @examples
#' df <- data.frame(
#'   Group = rep(c("R", "NR"), each = 30),
#'   CellType = rep(rep(c("DC", "CD8"), each = 15), 2),
#'   value = c(rnorm(15, 0.2, 0.05), rnorm(15, 0.35, 0.05),
#'             rnorm(15, 0.15, 0.05), rnorm(15, 0.3, 0.05))
#' )
#' plot_group_comparison(df, value_col = "value", group_col = "Group",
#'                       by = "CellType")
plot_groups <- function(data,
                                  value_col,
                                  group_col,
                                  by = NULL,
                                  method = "wilcox.test",
                                  comparisons = NULL,
                                  geom = c("both", "box", "violin"),
                                  percent = TRUE,
                                  y_lab = "% of DC 20\u00b5m neighbors",
                                  palette = NULL,
                                  jitter = FALSE,
                                  point_size = 0.8,
                                  alpha = 0.5,
                                  title = NULL,
                                  iqr_filter = TRUE) {
  vc <- .find_col(data, value_col)
  gc <- .find_col(data, group_col)
  geom <- match.arg(geom)
  if (is.null(vc) || is.null(gc)) {
    stop("Could not find the value and/or group column.", call. = FALSE)
  }
  bc <- if (is.null(by)) NULL else .find_col(data, by)

  d <- data
  if (isTRUE(iqr_filter)) {
    filter_group <- if (is.null(bc)) rep("all", nrow(d)) else d[[bc]]
    keep <- stats::ave(d[[vc]], filter_group, FUN = function(z) {
      q1 <- stats::quantile(z, 0.25, na.rm = TRUE)
      q3 <- stats::quantile(z, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      !is.na(z) & z >= q1 - 1.5 * iqr & z <= q3 + 1.5 * iqr
    })
    d <- d[as.logical(keep), , drop = FALSE]
  }
  if (percent) d[[vc]] <- d[[vc]] * 100

  groups <- unique(d[[gc]])
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2L) {
    stop("Need at least two groups in `group_col`.", call. = FALSE)
  }
  if (is.null(palette)) {
    if (length(groups) == 2L) {
      palette <- c(R = "#E64B35", NR = "#4DBBD5")
    } else {
      palette <- .default_palette(length(groups))
    }
  }
  if (is.null(names(palette))) {
    palette <- stats::setNames(palette, as.character(groups))
  }

  # colour (group outline) is only mapped when a violin is drawn; boxplots are
  # white so they only use the fill aesthetic
  use_color <- geom %in% c("both", "violin")
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[gc]], y = .data[[vc]], fill = .data[[gc]]
  ))
  if (use_color) {
    p <- p + ggplot2::aes(color = .data[[gc]])
  }
  if (geom %in% c("both", "violin")) {
    p <- p + ggplot2::geom_violin()
  }
  if (geom %in% c("both", "box")) {
    p <- p + ggplot2::geom_boxplot(
      width = 0.2, color = "white", alpha = 0.5, outlier.shape = NA
    )
  }
  if (!is.null(bc)) {
    p <- p + ggplot2::facet_wrap(
      stats::as.formula(paste("~", bc)), nrow = 1, scales = "free_y"
    )
  }
  p <- p +
    ggthemes::theme_few() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10, color = "black"),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      panel.border = ggplot2::element_rect(
        fill = NA, color = "black", linewidth = 1, linetype = "solid"
      ),
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 9, face = "bold")
    ) +
    ggplot2::labs(x = "", y = y_lab, title = title) +
    ggplot2::scale_fill_manual(values = palette)
  if (use_color) {
    p <- p + ggplot2::scale_color_manual(values = palette)
  }

  if (percent) {
    p <- p + ggplot2::scale_y_continuous(breaks = seq(0, 100, by = 2.5))
  }
  if (jitter) {
    p <- p + ggplot2::geom_jitter(width = 0.15, size = point_size, alpha = alpha)
  }
  if (is.null(comparisons)) {
    comparisons <- utils::combn(as.character(groups), 2, simplify = FALSE)
  }
  p_data <- d
  if (!is.null(bc)) {
    valid_levels <- unique(d[[bc]])[vapply(unique(d[[bc]]), function(level) {
      z <- d[d[[bc]] == level, , drop = FALSE]
      group_sizes <- table(z[[gc]])
      length(group_sizes) >= 2L && all(group_sizes >= 2L) &&
        length(unique(z[[vc]][!is.na(z[[vc]])])) > 1L
    }, logical(1))]
    p_data <- d[d[[bc]] %in% valid_levels, , drop = FALSE]
  }
  if (nrow(p_data) > 0L && length(unique(p_data[[vc]][!is.na(p_data[[vc]])])) > 1L) {
    p <- p + ggpubr::stat_compare_means(
      data = p_data,
      comparisons = comparisons,
      method = method,
      label = "p.format",
      size = 3.5
    )
  }
  p
}

#' @rdname plot_groups
#' @export
plot_group_comparison <- plot_groups
