#' Plot cell-type proportions per sample (stacked bars)
#'
#' Stacked bar chart showing, for each sample, either the proportion or the
#' absolute count of each cell type.
#'
#' @param data A data frame with a cell-type column and a sample column.
#' @param celltype_col Name of the cell-type column.
#' @param sample_col Name of the sample column.
#' @param levels Optional ordering of cell types (bottom to top of the bars).
#' @param proportion Logical; if `TRUE` (default) the y axis shows the
#'   proportion within each sample, otherwise absolute counts.
#' @param palette Optional character vector of colours.
#' @param title Optional plot title.
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
#'
#' @examples
#' data(mif_cells)
#' data(mif_cutoffs)
#' data(mif_rules)
#' marker_map <- c(
#'   CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
#'   CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
#'   PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
#' )
#' annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
#' plot_celltype_proportions(annotated)
plot_proportions <- function(data,
                                      celltype_col = "CellType",
                                      sample_col = "sample_id",
                                      levels = NULL,
                                      proportion = TRUE,
                                      palette = NULL,
                                      title = NULL) {
  cc <- .find_col(data, celltype_col)
  sc <- .find_col(data, sample_col)
  if (is.null(cc) || is.null(sc)) {
    stop("Could not find the cell-type and/or sample column.",
         call. = FALSE)
  }
  if (is.null(levels)) levels <- unique(data[[cc]])
  d <- data |>
    dplyr::mutate(.ct = factor(.data[[cc]], levels = levels)) |>
    dplyr::count(.data[[sc]], .ct, .drop = FALSE) |>
    dplyr::group_by(.data[[sc]]) |>
    dplyr::mutate(
      .total = sum(n),
      .prop = n / .total
    ) |>
    dplyr::ungroup()

  y_lab <- "Proportion" 
  if (proportion) {
    p <- ggplot2::ggplot(d, ggplot2::aes(
      x = .data[[sc]], y = .prop, fill = .ct
    ))
  } else {
    y_lab <- "Cell count"
    p <- ggplot2::ggplot(d, ggplot2::aes(
      x = .data[[sc]], y = n, fill = .ct
    ))
  }
  p <- p +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = y_lab, title = title)
  n_types <- length(unique(d$.ct))
  if (is.null(palette)) palette <- .default_palette(n_types)
  p <- p + ggplot2::scale_fill_manual(values = palette)
  if (proportion) {
    p <- p + ggplot2::scale_y_continuous(labels = scales::percent)
  }
  p
}

#' @rdname plot_proportions
#' @export
plot_celltype_proportions <- plot_proportions
