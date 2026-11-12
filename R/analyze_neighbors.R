#' Run the neighbour-analysis pipeline
#'
#' Convenience wrapper that runs [find_neighbors()], [neighbor_fov_summary()]
#' and [nearest_neighbor_rank()] in one call.
#'
#' @param data A data frame with FOV, cell-type and centroid columns.
#' @param radius Search radius in micrometres.
#' @param center_types Cell types used as centre cells for the FOV summary.
#' @inheritParams find_neighbors
#' @param sample_col Name of the sample column.
#'
#' @return A list with elements:
#' \itemize{
#'   \item `neighbors`: per-cell neighbour counts (output of
#'     [find_neighbors()]);
#'   \item `fov_summary`: per-FOV neighbour statistics (output of
#'     [neighbor_fov_summary()]);
#'   \item `nearest_rank`: per-FOV nearest-neighbour ranking (output of
#'     [nearest_neighbor_rank()]).
#' }
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
#' res <- analyze_neighbors(fov, radius = 20,
#'                          center_types = c("cDC", "cDC(Mature/Active)"))
#' names(res)
analyze_neighbors <- function(data,
                              radius = 20,
                              center_types = NULL,
                              celltype_col = "CellType",
                              fov_col = "FOV",
                              sample_col = "sample_id",
                              x = NULL,
                              y = NULL,
                              id_col = NULL,
                              max_k = 1000,
                              include_self = FALSE) {
  nb <- find_neighbors(
    data, radius = radius, celltype_col = celltype_col, fov_col = fov_col,
    x = x, y = y, id_col = id_col, max_k = max_k,
    include_self = include_self
  )
  summ <- neighbor_fov_summary(
    nb, center_types = center_types, celltype_col = celltype_col,
    fov_col = fov_col, sample_col = sample_col
  )
  rank <- nearest_neighbor_rank(
    nb, fov_col = fov_col, sample_col = sample_col,
    center_types = center_types, celltype_col = celltype_col
  )
  list(neighbors = nb, fov_summary = summ, nearest_rank = rank)
}
