#' Rank nearest-neighbour cell types per FOV
#'
#' For each FOV, tabulates the `nearestNGB` values (the cell type of the
#' closest neighbour) and ranks them by count (mirroring the nearest-neighbour
#' ranking of `mIF_Tzone_04_statistic.R`).
#'
#' @param data Output of [find_neighbors()].
#' @param fov_col Name of the FOV column.
#' @param sample_col Name of the sample column.
#' @param nearest_col Name of the nearest-neighbour column.
#' @param center_types Optional cell types to retain as centre cells. If
#'   `NULL`, all cells are used.
#' @param celltype_col Name of the centre-cell type column.
#' @param exclude Character vector of `nearestNGB` values to drop (default
#'   `"None"` for cells without neighbours).
#'
#' @return A data frame with columns FOV, Sample, `nearest_col`, `Count`,
#'   `Proportion` and `Rank` (1 = most frequent), sorted by FOV and rank.
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
#' nb <- find_neighbors(fov, radius = 60)
#' rank <- nearest_neighbor_rank(nb)
#' head(rank)
nearest_neighbor_rank <- function(data,
                                  fov_col = "FOV",
                                  sample_col = "sample_id",
                                  nearest_col = "nearestNGB",
                                  center_types = NULL,
                                  celltype_col = "CellType",
                                  exclude = "None") {
  fc <- .find_col(data, fov_col)
  sc <- .find_col(data, sample_col)
  nc <- .find_col(data, nearest_col)
  cc <- .find_col(data, celltype_col)
  if (is.null(fc) || is.null(sc) || is.null(nc) ||
      (!is.null(center_types) && is.null(cc))) {
    stop("Could not find the FOV, sample and/or nearest-neighbour columns.",
         call. = FALSE)
  }
  keep <- !is.na(data[[nc]]) & !data[[nc]] %in% exclude
  if (!is.null(center_types)) keep <- keep & data[[cc]] %in% center_types
  d <- data[keep, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No nearest-neighbour records after exclusion.", call. = FALSE)
  }
  d |>
    dplyr::mutate(.__sample = .data[[sc]]) |>
    dplyr::count(.data[[fc]], .__sample, .data[[nc]], name = "Count") |>
    dplyr::group_by(.data[[fc]]) |>
    dplyr::mutate(
      Proportion = Count / sum(Count),
      Rank = rank(-Count, ties.method = "min")
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data[[fc]], Rank) |>
    dplyr::rename(Sample = .__sample)
}
