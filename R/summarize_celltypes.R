#' Summarise cell-type counts per sample
#'
#' Builds a per-sample count table with one column per cell type and a
#' `Total` column. Cell-type columns are ordered by `levels` (or, by default,
#' by the order of appearance in `data`).
#'
#' @param data A data frame with a cell-type column and a sample column
#'   (e.g. the output of [run_annotation()]).
#' @param celltype_col Name of the cell-type column.
#' @param sample_col Name of the sample column.
#' @param levels Optional character vector giving the order of cell types in
#'   the output. If `NULL`, the order of appearance is used.
#' @param include_total Logical; if `TRUE` (default) add a `Total` column.
#'
#' @return A data frame with one row per sample.
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
#' summarize_celltypes(annotated)
summarize_celltypes <- function(data,
                                celltype_col = "CellType",
                                sample_col = "sample_id",
                                levels = NULL,
                                include_total = TRUE) {
  cc <- .find_col(data, celltype_col)
  sc <- .find_col(data, sample_col)
  if (is.null(cc)) {
    stop("Cell-type column '", celltype_col, "' not found in `data`.",
         call. = FALSE)
  }
  if (is.null(sc)) {
    stop("Sample column '", sample_col, "' not found in `data`.", call. = FALSE)
  }
  if (is.null(levels)) {
    levels <- unique(data[[cc]])
  } else {
    miss <- setdiff(levels, unique(data[[cc]]))
    if (length(miss)) {
      warning("Cell types in `levels` not present in `data`: ",
              paste(miss, collapse = ", "), call. = FALSE)
    }
    extra <- setdiff(unique(data[[cc]]), levels)
    if (length(extra)) levels <- c(levels, extra)
  }

  counts <- data |>
    dplyr::mutate(
      .ct_fct = factor(.data[[cc]], levels = levels),
      .sm_fct = as.character(.data[[sc]])
    ) |>
    dplyr::count(.sm_fct, .ct_fct, .drop = FALSE) |>
    tidyr::pivot_wider(names_from = .ct_fct, values_from = n,
                       values_fill = 0L) |>
    dplyr::rename(Sample = .sm_fct)

  if (include_total) {
    counts <- counts |>
      dplyr::mutate(Total = rowSums(
        dplyr::across(dplyr::all_of(levels)), na.rm = TRUE
      ))
  }
  counts
}