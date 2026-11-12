#' Keep only FOVs that pass quality control
#'
#' Filters cell-level (or summary-level) data down to the FOVs whose
#' `LowCellFlag` is one of `keep` (default `"OK"`).
#'
#' @param data Cell-level data (or a QC summary) with an FOV column.
#' @param qc The QC summary returned by [qc_fov_summary()].
#' @param keep Character vector of `LowCellFlag` values to keep.
#' @param fov_col Name of the FOV column in `data` and `qc`.
#'
#' @return The filtered data frame.
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
#' qc <- qc_fov_summary(fov)
#' clean <- filter_qc_fovs(fov, qc)
filter_qc_fovs <- function(data,
                           qc,
                           keep = "OK",
                           fov_col = "FOV") {
  fc <- .find_col(data, fov_col)
  qfc <- .find_col(qc, fov_col)
  flag <- .find_col(qc, "LowCellFlag")
  if (is.null(fc) || is.null(qfc)) {
    stop("Could not find the FOV column in `data` and/or `qc`.", call. = FALSE)
  }
  if (is.null(flag)) {
    stop("`qc` must contain a LowCellFlag column (see qc_fov_summary()).",
         call. = FALSE)
  }
  keep_fovs <- qc[[qfc]][qc[[flag]] %in% keep]
  data[data[[fc]] %in% keep_fovs, , drop = FALSE]
}