#' Run the full Phase-1 annotation pipeline
#'
#' Convenience wrapper that chains [filter_region()] (optional),
#' [apply_cutoffs()] and [annotate_cells()] on a single data frame.
#'
#' @param data A data frame with QuPath cell data (must contain a sample
#'   column and the intensity columns named by `marker_map`).
#' @param cutoffs A cutoff data frame (see [apply_cutoffs()]).
#' @param rules A rule data frame (see [annotate_cells()]).
#' @param marker_map Named character vector mapping markers to intensity
#'   columns (see [apply_cutoffs()]).
#' @param regions Optional character vector of parent regions to keep. If
#'   `NULL`, no region filtering is performed.
#' @param parent_col Name of the parent column (see [filter_region()]).
#' @param region_match Region matching mode passed to [filter_region()].
#' @param region_ignore_case Logical; ignore case when matching region names.
#' @param sample_col Name of the sample column in `data`.
#' @param celltype_col Name of the output cell-type column.
#' @param unmatched Label for cells matching no rule (see [annotate_cells()]).
#' @param ... Further arguments passed to [annotate_cells()].
#'
#' @return The annotated data frame.
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
run_annotation <- function(data,
                           cutoffs,
                           rules,
                           marker_map,
                           regions = NULL,
                           parent_col = NULL,
                           region_match = c("exact", "fixed", "regex"),
                           region_ignore_case = FALSE,
                           sample_col = "sample_id",
                           celltype_col = "CellType",
                           unmatched = "Other",
                           ...) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  d <- data
  if (!is.null(regions)) {
    d <- filter_region(
      d, parent_col = parent_col, regions = regions,
      match = match.arg(region_match), ignore_case = region_ignore_case
    )
  }
  d <- apply_cutoffs(d, cutoffs, marker_map, sample_col = sample_col)
  d <- annotate_cells(d, rules,
                      celltype_col = celltype_col,
                      unmatched = unmatched,
                      ...)
  d
}
