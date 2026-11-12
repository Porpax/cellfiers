#' Save one cell map per sample
#'
#' Loops over samples and saves an individual [plot_cell_map()] figure per
#' sample into `output_dir`, named `<sample>_cell_map.<filetype>`.
#'
#' @param data A data frame with sample, centroid and colour columns.
#' @param output_dir Directory where figures are written (created if needed).
#' @param filetype Output file type passed to [ggplot2::ggsave()] (e.g.
#'   `"png"` or `"pdf"`).
#' @param width,height,dpi Figure size in inches and resolution.
#' @param ... Further arguments passed to [plot_cell_map()] (e.g. `color`,
#'   `palette`, `point_size`).
#'
#' @return Invisibly, the file paths written.
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
#' \dontrun{
#' plot_cell_maps_by_sample(annotated, "figures")
#' }
plot_cells_by_sample <- function(data,
                                     output_dir,
                                     filetype = "png",
                                     width = 8,
                                     height = 6,
                                     dpi = 150,
                                     ...) {
  sc <- .find_col(data, "sample_id")
  if (is.null(sc)) {
    stop("Could not find the sample column. Pass data with a sample column.",
         call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  samples <- unique(data[[sc]])
  paths <- character(length(samples))
  for (i in seq_along(samples)) {
    s <- samples[i]
    p <- plot_cell_map(
      data[data[[sc]] == s, , drop = FALSE],
      sample_id = s,
      ...
    )
    paths[i] <- file.path(output_dir, paste0(s, "_cell_map.", filetype))
    ggplot2::ggsave(paths[i], p, width = width, height = height, dpi = dpi)
  }
  invisible(paths)
}

#' Save one FOV distribution plot per sample
#'
#' Loops over samples and saves an individual [plot_fov_distribution()]
#' figure per sample into `output_dir`, named
#' `<sample>_FOV_distribution.<filetype>`.
#'
#' @param data A data frame with sample, FOV and centroid columns (output of
#'   [assign_fov()]).
#' @param output_dir Directory where figures are written (created if needed).
#' @param filetype Output file type passed to [ggplot2::ggsave()].
#' @param width,height,dpi Figure size in inches and resolution.
#' @param ... Further arguments passed to [plot_fov_distribution()].
#'
#' @return Invisibly, the file paths written.
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
#' fov <- assign_fov(annotated)
#' \dontrun{
#' plot_fov_distributions_by_sample(fov, "figures")
#' }
plot_fovs_by_sample <- function(data,
                                             output_dir,
                                             filetype = "png",
                                             width = 8,
                                             height = 6,
                                             dpi = 150,
                                             ...) {
  sc <- .find_col(data, "sample_id")
  if (is.null(sc)) {
    stop("Could not find the sample column. Pass data with a sample column.",
         call. = FALSE)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  samples <- unique(data[[sc]])
  paths <- character(length(samples))
  for (i in seq_along(samples)) {
    s <- samples[i]
    p <- plot_fov_distribution(
      data[data[[sc]] == s, , drop = FALSE],
      sample_id = s,
      ...
    )
    paths[i] <- file.path(output_dir,
                          paste0(s, "_FOV_distribution.", filetype))
    ggplot2::ggsave(paths[i], p, width = width, height = height, dpi = dpi)
  }
  invisible(paths)
}

#' @rdname plot_cells_by_sample
#' @export
plot_cell_maps_by_sample <- plot_cells_by_sample

#' @rdname plot_fovs_by_sample
#' @export
plot_fov_distributions_by_sample <- plot_fovs_by_sample
