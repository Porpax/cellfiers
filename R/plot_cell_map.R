# encoding: UTF-8
#' Plot a spatial map of cells coloured by cell type
#'
#' Scatter plot of cell centroids, coloured by a categorical column (cell
#' type by default). With multiple samples the panels are faceted by sample,
#' or pass `sample_id` to plot a single sample. This is the "display"
#' companion of the annotation step and helps to visually inspect where each
#' cell type sits in the tissue.
#'
#' @param data A data frame with x/y coordinates and a colour column.
#' @param x,y Names of the centroid columns. Defaults to QuPath-style
#'   `"Centroid X µm"` / `"Centroid Y µm"` (or normalised matches).
#' @param color Name of the column used to colour the points. Defaults to a
#'   column named `"CellType"`, then `"Name"`, then `"Classification"`.
#' @param sample_col Name of the sample column used for faceting/filtering.
#' @param sample_id Optional sample id; if given, only that sample is plotted
#'   (and used as the default title).
#' @param facet Logical; if `TRUE` (default) facet by `sample_col` when more
#'   than one sample is present and `sample_id` is `NULL`.
#' @param point_size Point size passed to [ggplot2::geom_point()].
#' @param alpha Point transparency passed to [ggplot2::geom_point()].
#' @param palette Optional character vector of colours; by default a
#'   colour-blind friendly palette is used.
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
#' plot_cell_map(annotated)
#' plot_cell_map(annotated, sample_id = "C1001")
plot_cell_map <- function(data,
                          x = NULL,
                          y = NULL,
                          color = NULL,
                          sample_col = "sample_id",
                          sample_id = NULL,
                          facet = TRUE,
                          point_size = 0.6,
                          alpha = 0.6,
                          palette = NULL,
                          title = NULL) {
  xc <- .find_col(data, c(x, "Centroid X um", "Centroid X", "x", "X"))
  yc <- .find_col(data, c(y, "Centroid Y um", "Centroid Y", "y", "Y"))
  cc <- .find_col(data, c(color, "CellType", "Name", "Classification"))
  sc <- .find_col(data, sample_col)
  if (is.null(xc) || is.null(yc)) {
    stop("Could not find x/y centroid columns. Pass `x` and `y` explicitly.",
         call. = FALSE)
  }

  if (!is.null(sample_id)) {
    if (is.null(sc)) {
      stop("`sample_id` was given but no sample column was found.",
           call. = FALSE)
    }
    data <- data[data[[sc]] == sample_id, , drop = FALSE]
    if (nrow(data) == 0L) {
      stop("No cells found for sample_id '", sample_id, "'.", call. = FALSE)
    }
    if (is.null(title)) title <- sample_id
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[[xc]], y = .data[[yc]]
  ))
  if (is.null(cc)) {
    p <- p + ggplot2::geom_point(size = point_size, alpha = alpha)
  } else {
    p <- p + ggplot2::aes(color = .data[[cc]]) +
      ggplot2::geom_point(size = point_size, alpha = alpha)
    n_types <- length(unique(data[[cc]]))
    if (is.null(palette)) palette <- .default_palette(n_types)
    p <- p + ggplot2::scale_color_manual(values = palette)
  }
  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(
      x = xc, y = yc,
      color = if (is.null(cc)) NULL else cc,
      title = title
    )
  if (facet && is.null(sample_id) && !is.null(sc) &&
      dplyr::n_distinct(data[[sc]]) > 1) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", sc)))
  }
  p
}