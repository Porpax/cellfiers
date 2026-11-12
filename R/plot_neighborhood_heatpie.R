# encoding: UTF-8
#' Plot an FOV heat-pie map for centre--target neighbourhoods
#'
#' Tile colour represents mean target-type richness, pie sectors represent
#' weighted target composition, and pie area represents mean weighted target
#' abundance per centre. By default, target-rich FOVs use the project's
#' pale-pink to deep-rose palette.
#'
#' @param result Output from [target_neighbors()].
#' @param target_colors Named colours for target types. Names may be target
#'   type labels or the corresponding `Weighted_...` columns.
#' @param richness_colors Colours from low to high richness.
#' @param other_selected_color Colour for FOVs dominated by user-specified
#'   non-target types.
#' @param centre_only_color Colour for centre--centre FOVs without targets.
#' @param other_isolated_color Colour for other/isolated centre FOVs.
#' @param no_centre_color Colour for FOVs without a centre cell.
#' @param border_color Grid-border colour.
#' @param max_richness Upper limit of the richness legend. Defaults to the
#'   number of requested target types.
#' @param pie_max Optional common upper reference for pie size.
#' @param pie_radius_fraction Maximum pie radius as a fraction of grid size.
#' @param title,subtitle Optional plot text.
#'
#' @return A ggplot object.
#' @examples
#' cells <- data.frame(
#'   sample_id = "S1", FOV = "F1", CellType = c("DC", "CD8 T"),
#'   x = c(0, 10), y = c(0, 0)
#' )
#' result <- target_neighbors(cells, "DC", "CD8 T", radius = 20)
#' plot_heatpie(result)
#' @export
plot_heatpie <- function(
    result,
    target_colors = NULL,
    richness_colors = c("#FECEF5", "#FECEF5", "#EC407A", "#B6066F"),
    other_selected_color = "#7ABD31",
    centre_only_color = "#FFFFCC",
    other_isolated_color = "#CBCBCB",
    no_centre_color = "white",
    border_color = "#666666",
    max_richness = NULL,
    pie_max = NULL,
    pie_radius_fraction = 0.38,
    title = NULL,
    subtitle = NULL) {
  if (!is.list(result) || is.null(result$fov_level) ||
      is.null(result$weighted_columns) || is.null(result$settings)) {
    stop("`result` must be returned by `target_neighbors()`.",
         call. = FALSE)
  }
  fov <- result$fov_level
  weighted_columns <- result$weighted_columns
  missing_columns <- setdiff(
    c("FOV_center_X", "FOV_center_Y", "FOVClass", "MeanTargetRichness",
      "MeanWeightedTargetAbundancePerCenter", weighted_columns), names(fov)
  )
  if (length(missing_columns)) {
    stop("The result is missing required columns: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  grid_size <- result$settings$grid_size
  target_types <- result$settings$target_types
  if (is.null(max_richness)) max_richness <- length(target_types)
  if (!is.numeric(max_richness) || max_richness <= 0) {
    stop("`max_richness` must be positive.", call. = FALSE)
  }

  if (is.null(target_colors)) {
    target_colors <- stats::setNames(
      scales::hue_pal()(length(target_types)), weighted_columns
    )
  } else {
    if (is.null(names(target_colors))) {
      if (length(target_colors) != length(target_types)) {
        stop("Unnamed `target_colors` must match `target_types` in length.",
             call. = FALSE)
      }
      names(target_colors) <- weighted_columns
    } else if (all(target_types %in% names(target_colors))) {
      target_colors <- stats::setNames(
        unname(target_colors[target_types]), weighted_columns
      )
    } else {
      target_colors <- target_colors[weighted_columns]
    }
    if (anyNA(target_colors)) {
      stop("`target_colors` must provide a colour for every target type.",
           call. = FALSE)
    }
  }

  pie_data <- fov[
    fov$N_Center > 0 & fov$MeanWeightedTargetAbundancePerCenter > 0,
    , drop = FALSE
  ]
  if (is.null(pie_max)) {
    pie_max <- if (nrow(pie_data)) {
      max(pie_data$MeanWeightedTargetAbundancePerCenter, na.rm = TRUE)
    } else {
      1
    }
  }
  if (!is.finite(pie_max) || pie_max <= 0) pie_max <- 1
  pie_data$PieRadius_um <- pie_radius_fraction * grid_size * sqrt(
    pmin(pie_data$MeanWeightedTargetAbundancePerCenter / pie_max, 1)
  )

  p <- ggplot2::ggplot(
    fov, ggplot2::aes(x = .data$FOV_center_X, y = .data$FOV_center_Y)
  ) +
    ggplot2::geom_tile(
      data = fov[fov$FOVClass == "No-centre", , drop = FALSE],
      width = grid_size, height = grid_size,
      fill = no_centre_color, color = border_color, linewidth = 0.42
    ) +
    ggplot2::geom_tile(
      data = fov[fov$FOVClass == "Other-selected", , drop = FALSE],
      width = grid_size, height = grid_size,
      fill = other_selected_color, color = border_color, linewidth = 0.42
    ) +
    ggplot2::geom_tile(
      data = fov[fov$FOVClass == "Centre-only", , drop = FALSE],
      width = grid_size, height = grid_size,
      fill = centre_only_color, color = border_color, linewidth = 0.42
    ) +
    ggplot2::geom_tile(
      data = fov[fov$FOVClass == "Other/isolated", , drop = FALSE],
      width = grid_size, height = grid_size,
      fill = other_isolated_color, color = border_color, linewidth = 0.42
    ) +
    ggplot2::geom_tile(
      data = fov[fov$FOVClass == "Target", , drop = FALSE],
      ggplot2::aes(fill = .data$MeanTargetRichness),
      width = grid_size, height = grid_size,
      color = border_color, linewidth = 0.42
    ) +
    ggplot2::scale_fill_gradientn(
      colors = richness_colors,
      values = c(0, min(1 / max_richness, 0.4), 0.65, 1),
      limits = c(0, max_richness), breaks = 0:max_richness,
      oob = scales::squish, na.value = no_centre_color,
      name = "Mean target richness"
    )

  if (nrow(pie_data)) {
    p <- p +
      ggnewscale::new_scale_fill() +
      scatterpie::geom_scatterpie(
        data = pie_data,
        ggplot2::aes(
          x = .data$FOV_center_X, y = .data$FOV_center_Y,
          r = .data$PieRadius_um
        ),
        cols = weighted_columns,
        color = "#303030", linewidth = 0.28, alpha = 0.88
      ) +
      ggplot2::scale_fill_manual(
        values = target_colors,
        labels = stats::setNames(target_types, weighted_columns),
        name = "Weighted target composition"
      )
  }

  if (is.null(title)) {
    title <- paste0(
      result$settings$sample, ": ",
      paste(result$settings$center_types, collapse = " + "),
      "-centred ", result$settings$radius, "-um neighbourhoods"
    )
  }
  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Pink-rose = target richness; pale yellow = centre-centre; ",
      "pie area = mean weighted target abundance/centre; weighting = ",
      result$settings$shared_weight
    )
  }
  p +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill = NA, color = "#444444", linewidth = 0.7
      ),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Centroid X (um)", y = "Centroid Y (um)"
    )
}

#' @rdname plot_heatpie
#' @export
plot_neighborhood_heatpie <- plot_heatpie
