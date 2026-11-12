#' Annotate a directory of QuPath exports and save results
#'
#' Batch pipeline over a directory of QuPath classification exports: reads
#' each file, optionally filters to a region, applies cutoffs and rules,
#' writes one annotated CSV per sample plus a combined cell-type summary to
#' `output_dir`.
#'
#' @param input_dir Directory containing QuPath export files.
#' @param output_dir Directory where annotated CSVs and the summary are
#'   written (created if needed).
#' @param cutoffs A cutoff data frame (see [apply_cutoffs()]).
#' @param rules A rule data frame (see [annotate_cells()]).
#' @param marker_map Named character vector mapping markers to intensity
#'   columns (see [apply_cutoffs()]).
#' @param regions Optional character vector of parent regions to keep.
#' @param pattern File-name pattern passed to [read_qupath_batch()].
#' @param sample_pattern Regex with one capture group to extract sample ids
#'   (see [read_qupath()]).
#' @param encoding File encoding passed to [readr::read_tsv()].
#' @param celltype_col Name of the output cell-type column.
#' @param write_summary Logical; if `TRUE` (default) also write
#'   `celltype_summary.csv` with per-sample cell-type counts.
#' @param ... Further arguments passed to [annotate_cells()].
#'
#' @return Invisibly, a named list of annotated data frames (one per sample).
#' @export
#'
#' @examples
#' \dontrun{
#' annotate_files(
#'   input_dir = "classfy_data",
#'   output_dir = "annotated_results",
#'   cutoffs = cutoffs_df,
#'   rules = rules_df,
#'   marker_map = marker_map,
#'   regions = c("Annotation (SLN_Tregion)", "SLN_Tregion")
#' )
#' }
annotate_files <- function(input_dir,
                           output_dir,
                           cutoffs,
                           rules,
                           marker_map,
                           regions = NULL,
                           pattern = "[.]txt$",
                           sample_pattern = ".*_(C[0-9]+)[.]txt$",
                           encoding = "UTF-8",
                           celltype_col = "CellType",
                           write_summary = TRUE,
                           ...) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_list <- list()
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No files matching pattern ", pattern, " found in ", input_dir,
         call. = FALSE)
  }
  for (f in files) {
    d <- read_qupath(f, sample_pattern = sample_pattern, encoding = encoding)
    d <- run_annotation(
      d, cutoffs, rules, marker_map,
      regions = regions,
      celltype_col = celltype_col,
      ...
    )
    sid <- d[["sample_id"]][1]
    out_file <- file.path(output_dir, paste0(sid, "_annotated.csv"))
    readr::write_csv(d, out_file)
    message("Annotated ", sid, ": ", nrow(d), " cells -> ", out_file)
    out_list[[sid]] <- d
  }
  if (write_summary && length(out_list)) {
    combined <- dplyr::bind_rows(out_list)
    summ <- summarize_celltypes(
      combined,
      celltype_col = celltype_col,
      sample_col = "sample_id"
    )
    summ_file <- file.path(output_dir, "celltype_summary.csv")
    utils::write.csv(summ, summ_file, row.names = FALSE)
    message("Wrote summary -> ", summ_file)
  }
  invisible(out_list)
}