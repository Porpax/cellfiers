#' Plot grouped neighbour proportions as a ring radar chart
#'
#' Draws one closed radar line per group and an outer segmented annotation
#' ring. The function is designed for sample-level output from
#' [neighbor_sample_summary()], but accepts any long data frame containing a
#' group, neighbour cell type, and numeric proportion.
#'
#' @param data Long-format data frame.
#' @param group_col Name of the grouping column.
#' @param category_col Name of the neighbour cell-type/category column.
#' @param value_col Name of the numeric proportion column.
#' @param group_colors Colours for group lines and points. A named vector is
#'   recommended, for example `c(R = "#C44732", NR = "#087DB5")`.
#' @param ring_colors One outer-ring colour, an unnamed vector in category
#'   order, or a named vector keyed by category.
#' @param category_order Optional category order around the circle.
#' @param category_labels Optional named character vector replacing displayed
#'   category labels without changing the data.
#' @param p_values Optional data frame containing category and p-value columns.
#'   When omitted and `show_p = TRUE`, two-group Wilcoxon tests are calculated
#'   from the rows in `data`.
#' @param p_category_col,p_value_col Column names in `p_values`.
#' @param show_p Logical; display one p-value per category.
#' @param breaks Numeric radial grid breaks. If `NULL`, pretty breaks are used.
#' @param limits Numeric radial limits for the data region. If `NULL`, limits
#'   are inferred from zero and the observed values.
#' @param ring_width Width of the outer ring relative to the data radius.
#' @param ring_alpha Outer-ring opacity.
#' @param line_width,point_size Radar line width and point size.
#' @param label_size,p_size,grid_label_size Text sizes for ring labels,
#'   p-values, and radial tick labels.
#' @param title Optional plot title.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#'
#' @examples
#' d <- expand.grid(
#'   Sample = paste0("S", 1:6),
#'   CellType = c("T cell", "DC", "Other")
#' )
#' d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
#' d$mean_prop <- seq(0.1, 0.6, length.out = nrow(d))
#' plot_neighbor_ring_radar(
#'   d, group_col = "Group", category_col = "CellType",
#'   value_col = "mean_prop",
#'   group_colors = c(R = "#C44732", NR = "#087DB5"),
#'   ring_colors = "#D8C0EA"
#' )
plot_ring_radar <- function(
    data,
    group_col = "Group",
    category_col = "CellType",
    value_col = "mean_prop",
    group_colors = c(R = "#C44732", NR = "#087DB5"),
    ring_colors = "#D8C0EA",
    category_order = NULL,
    category_labels = NULL,
    p_values = NULL,
    p_category_col = category_col,
    p_value_col = "p.value",
    show_p = TRUE,
    breaks = NULL,
    limits = NULL,
    ring_width = 0.28,
    ring_alpha = 1,
    line_width = 1,
    point_size = 2.2,
    label_size = 3.6,
    p_size = 3,
    grid_label_size = 2.8,
    title = NULL) {
  gc <- .find_col(data, group_col)
  cc <- .find_col(data, category_col)
  vc <- .find_col(data, value_col)
  if (is.null(gc) || is.null(cc) || is.null(vc)) {
    stop("Could not find the group, category, and/or value column.",
         call. = FALSE)
  }
  if (!is.numeric(data[[vc]]) || !nrow(data)) {
    stop("`value_col` must be numeric and `data` must not be empty.",
         call. = FALSE)
  }
  d <- data[!is.na(data[[gc]]) & !is.na(data[[cc]]) &
              is.finite(data[[vc]]), , drop = FALSE]
  groups <- unique(as.character(d[[gc]]))
  categories <- if (is.null(category_order)) {
    unique(as.character(d[[cc]]))
  } else {
    as.character(category_order)
  }
  if (length(groups) < 1L || length(categories) < 3L) {
    stop("At least one group and three categories are required.",
         call. = FALSE)
  }
  missing_categories <- setdiff(categories, unique(as.character(d[[cc]])))
  if (length(missing_categories)) {
    stop("Categories absent from `data`: ",
         paste(missing_categories, collapse = ", "), call. = FALSE)
  }

  radar <- d |>
    dplyr::mutate(
      .__group = as.character(.data[[gc]]),
      .__category = as.character(.data[[cc]])
    ) |>
    dplyr::filter(.data$.__category %in% categories) |>
    dplyr::group_by(.data$.__group, .data$.__category) |>
    dplyr::summarise(value = mean(.data[[vc]], na.rm = TRUE), .groups = "drop")
  complete <- expand.grid(
    .__group = groups, .__category = categories,
    stringsAsFactors = FALSE
  )
  radar <- dplyr::left_join(complete, radar,
                            by = c(".__group", ".__category"))
  if (anyNA(radar$value)) {
    stop("Every group must contain every requested category.", call. = FALSE)
  }
  radar$axis_id <- match(radar$.__category, categories)

  if (is.null(limits)) {
    limits <- c(0, max(radar$value, na.rm = TRUE))
    if (limits[2] <= 0) limits[2] <- 1
  }
  if (!is.numeric(limits) || length(limits) != 2L ||
      any(!is.finite(limits)) || limits[1] >= limits[2]) {
    stop("`limits` must contain two increasing finite numbers.",
         call. = FALSE)
  }
  if (any(radar$value < limits[1] | radar$value > limits[2])) {
    stop("Observed values fall outside `limits`.", call. = FALSE)
  }
  data_radius <- diff(limits)
  if (is.null(breaks)) {
    breaks <- pretty(limits, n = 5)
    breaks <- breaks[breaks >= limits[1] & breaks <= limits[2]]
  }
  ring_inner <- 1.13
  ring_outer <- ring_inner + ring_width
  label_radius <- (ring_inner + ring_outer) / 2
  p_radius <- 1.055
  axis_theta <- pi / 2 - 2 * pi * (seq_along(categories) - 1) /
    length(categories)
  radar$radius <- (radar$value - limits[1]) / data_radius
  radar$theta <- axis_theta[radar$axis_id]
  radar$plot_x <- radar$radius * cos(radar$theta)
  radar$plot_y <- radar$radius * sin(radar$theta)

  if (length(group_colors) < length(groups)) {
    stop("Supply at least one colour per group.", call. = FALSE)
  }
  if (is.null(names(group_colors))) {
    group_colors <- stats::setNames(group_colors[seq_along(groups)], groups)
  } else if (!all(groups %in% names(group_colors))) {
    stop("Named `group_colors` must include every group.", call. = FALSE)
  }
  if (length(ring_colors) == 1L) ring_colors <- rep(ring_colors, length(categories))
  if (!is.null(names(ring_colors))) {
    if (!all(categories %in% names(ring_colors))) {
      stop("Named `ring_colors` must include every category.", call. = FALSE)
    }
    ring_colors <- unname(ring_colors[categories])
  } else if (length(ring_colors) != length(categories)) {
    stop("`ring_colors` must have length 1 or one colour per category.",
         call. = FALSE)
  }

  ring <- data.frame(
    axis_id = seq_along(categories),
    category = categories,
    fill_key = categories,
    stringsAsFactors = FALSE
  )
  ring$label <- categories
  if (!is.null(category_labels)) {
    if (is.null(names(category_labels))) {
      stop("`category_labels` must be a named character vector.",
           call. = FALSE)
    }
    replace <- category_labels[ring$category]
    ring$label[!is.na(replace)] <- unname(replace[!is.na(replace)])
  }
  ring$theta <- axis_theta
  ring$label_x <- label_radius * cos(ring$theta)
  ring$label_y <- label_radius * sin(ring$theta)
  ring$angle <- ring$theta * 180 / pi - 90
  ring$angle[ring$angle > 90] <- ring$angle[ring$angle > 90] - 180
  ring$angle[ring$angle < -90] <- ring$angle[ring$angle < -90] + 180

  half_sector <- pi / length(categories)
  ring_polygons <- dplyr::bind_rows(lapply(seq_along(categories), function(i) {
    theta_outer <- seq(
      axis_theta[i] + half_sector,
      axis_theta[i] - half_sector,
      length.out = 60
    )
    theta_inner <- rev(theta_outer)
    data.frame(
      plot_x = c(ring_outer * cos(theta_outer),
                 ring_inner * cos(theta_inner)),
      plot_y = c(ring_outer * sin(theta_outer),
                 ring_inner * sin(theta_inner)),
      fill_key = categories[i],
      polygon_id = i
    )
  }))

  circle_theta <- seq(0, 2 * pi, length.out = 241)
  break_radius <- (breaks - limits[1]) / data_radius
  grid_circles <- dplyr::bind_rows(lapply(seq_along(break_radius), function(i) {
    data.frame(
      plot_x = break_radius[i] * cos(circle_theta),
      plot_y = break_radius[i] * sin(circle_theta),
      grid_id = i
    )
  }))
  spokes <- data.frame(
    x = 0, y = 0,
    xend = cos(axis_theta), yend = sin(axis_theta)
  )
  tick_theta <- pi / 2 - 0.035
  grid_labels <- data.frame(
    plot_x = break_radius * cos(tick_theta),
    plot_y = break_radius * sin(tick_theta),
    label = format(breaks, trim = TRUE, digits = 3)
  )

  p_labels <- NULL
  if (isTRUE(show_p)) {
    if (is.null(p_values)) {
      if (length(groups) != 2L) {
        stop("Automatic p-values require exactly two groups.", call. = FALSE)
      }
      p_num <- vapply(categories, function(category) {
        z <- d[as.character(d[[cc]]) == category, , drop = FALSE]
        x <- z[[vc]][as.character(z[[gc]]) == groups[1]]
        y <- z[[vc]][as.character(z[[gc]]) == groups[2]]
        x <- x[is.finite(x)]; y <- y[is.finite(y)]
        if (length(x) < 2L || length(y) < 2L ||
            length(unique(c(x, y))) < 2L) return(NA_real_)
        suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
      }, numeric(1))
      p_labels <- data.frame(category = categories, p.value = p_num)
    } else {
      pc <- .find_col(p_values, p_category_col)
      pv <- .find_col(p_values, p_value_col)
      if (is.null(pc) || is.null(pv) || !is.numeric(p_values[[pv]])) {
        stop("Could not find valid category and p-value columns in `p_values`.",
             call. = FALSE)
      }
      p_labels <- data.frame(
        category = as.character(p_values[[pc]]),
        p.value = p_values[[pv]], stringsAsFactors = FALSE
      )
      p_labels <- p_labels[match(categories, p_labels$category), , drop = FALSE]
    }
    p_labels$theta <- axis_theta
    p_labels$plot_x <- p_radius * cos(p_labels$theta)
    p_labels$plot_y <- p_radius * sin(p_labels$theta)
    p_labels$label <- ifelse(
      is.na(p_labels$p.value), "P = NA",
      paste0("P = ", format.pval(p_labels$p.value, digits = 2, eps = 2.2e-16))
    )
    p_labels$angle <- ring$angle
  }

  ring_scale <- stats::setNames(ring_colors, categories)
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = ring_polygons,
      ggplot2::aes(
        x = .data$plot_x, y = .data$plot_y,
        group = .data$polygon_id, fill = .data$fill_key
      ),
      color = NA, alpha = ring_alpha
    ) +
    ggplot2::geom_path(
      data = grid_circles,
      ggplot2::aes(
        x = .data$plot_x, y = .data$plot_y, group = .data$grid_id
      ),
      color = "grey82", linewidth = 0.65
    ) +
    ggplot2::geom_segment(
      data = spokes,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey82", linewidth = 0.65
    ) +
    ggplot2::geom_polygon(
      data = radar,
      ggplot2::aes(
        x = .data$plot_x, y = .data$plot_y,
        group = .data$.__group, color = .data$.__group
      ),
      fill = NA, linewidth = line_width
    ) +
    ggplot2::geom_point(
      data = radar,
      ggplot2::aes(
        x = .data$plot_x, y = .data$plot_y, color = .data$.__group
      ),
      size = point_size
    ) +
    ggplot2::geom_text(
      data = grid_labels,
      ggplot2::aes(
        x = .data$plot_x, y = .data$plot_y, label = .data$label
      ),
      size = grid_label_size, color = "grey25", hjust = -0.15,
      vjust = -0.15
    ) +
    ggplot2::geom_text(
      data = ring,
      ggplot2::aes(
        x = .data$label_x, y = .data$label_y,
        label = .data$label, angle = .data$angle
      ),
      fontface = "bold", size = label_size
    ) +
    ggplot2::scale_fill_manual(values = ring_scale, guide = "none") +
    ggplot2::scale_color_manual(values = group_colors, name = NULL) +
    ggplot2::coord_equal(
      xlim = c(-ring_outer * 1.06, ring_outer * 1.06),
      ylim = c(-ring_outer * 1.06, ring_outer * 1.06),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    ) +
    ggplot2::labs(title = title, color = NULL)

  if (!is.null(p_labels)) {
    p <- p + ggplot2::geom_text(
      data = p_labels,
      ggplot2::aes(
        x = .data$plot_x, y = .data$plot_y,
        label = .data$label, angle = .data$angle
      ),
      size = p_size
    )
  }
  p
}

#' @rdname plot_ring_radar
#' @export
plot_neighbor_ring_radar <- plot_ring_radar
