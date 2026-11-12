# encoding: UTF-8
#' Read a QuPath cell classification export
#'
#' Reads a single tab-separated classification export produced by QuPath
#' (e.g. Shanxi_SLN_C2309160.txt), adds a sample_id column and a unique
#' cell_id column (sample_id + "_" + row number). Column names are preserved
#' as they appear in the QuPath file (for example "Centroid X µm" or
#' "Cell: Opal 780 mean").
#'
#' @param file Path to a QuPath export file (tab-separated text).
#' @param sample_id Optional sample id. If NULL, it is extracted from the
#'   file name using sample_pattern.
#' @param sample_pattern Regular expression with exactly one capture group
#'   used to extract the sample id from the file name. Defaults to
#'   ".*_(C[0-9]+)[.]txt$" which matches names such as
#'   Shanxi_SLN_C2309160.txt and yields C2309160.
#' @param encoding File encoding passed to readr::read_tsv(). Defaults to
#'   "UTF-8".
#' @param ... Further arguments passed to readr::read_tsv().
#'
#' @return A tibble with the original QuPath columns plus sample_id and
#'   cell_id.
#' @export
#'
#' @examples
#' \dontrun{
#' cells <- read_qupath("classfy_data/Shanxi_SLN_C2309160.txt")
#' }
read_qupath <- function(file,
                        sample_id = NULL,
                        sample_pattern = ".*_(C[0-9]+)[.]txt$",
                        encoding = "UTF-8",
                        ...) {
  if (length(file) != 1L || !file.exists(file)) {
    stop("File does not exist: ", file, call. = FALSE)
  }
  if (is.null(sample_id)) {
    sample_id <- .extract_sample_id(file, sample_pattern)
  }
  df <- readr::read_tsv(
    file,
    locale = readr::locale(encoding = encoding),
    col_names = TRUE,
    trim_ws = TRUE,
    show_col_types = FALSE,
    ...
  )
  df[["sample_id"]] <- sample_id
  df <- .add_cell_id(df, "sample_id", "cell_id")
  df
}

#' Read all QuPath classification exports in a directory
#'
#' Batch version of read_qupath(). Reads every file matching pattern in dir,
#' extracts a sample id from each file name and returns either a list of data
#' frames (one per sample) or a single combined data frame.
#'
#' @param dir Path to a directory containing QuPath export files.
#' @param pattern File-name pattern passed to base::list.files().
#' @param sample_pattern Regular expression with one capture group used to
#'   extract sample ids (see read_qupath()).
#' @param combine Logical; if TRUE (default) all samples are combined into a
#'   single data frame with a sample_id column. If FALSE, a named list is
#'   returned instead.
#' @param encoding File encoding passed to readr::read_tsv().
#' @param skip_invalid Logical; if `TRUE` (default), files that cannot be
#'   read or whose columns have invalid (oversized/duplicated) names are
#'   skipped with a warning instead of aborting the batch.
#' @param ... Further arguments passed to read_qupath().
#'
#' @return A data frame (if combine = TRUE) or a named list of data frames.
#' @export
#'
#' @examples
#' \dontrun{
#' cells <- read_qupath_batch("classfy_data/")
#' }
read_qupath_batch <- function(dir,
                              pattern = "[.]txt$",
                              sample_pattern = ".*_(C[0-9]+)[.]txt$",
                              combine = TRUE,
                              encoding = "UTF-8",
                              skip_invalid = TRUE,
                              ...) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No files matching pattern ", pattern, " found in ", dir, call. = FALSE)
  }
  out <- list()
  skipped <- character()
  for (f in files) {
    d <- tryCatch(
      read_qupath(f, sample_pattern = sample_pattern, encoding = encoding, ...),
      error = function(e) e
    )
    if (inherits(d, "error")) {
      if (!skip_invalid) stop(d)
      warning("Skipping ", basename(f), ": ", conditionMessage(d), call. = FALSE)
      skipped <- c(skipped, basename(f))
      next
    }
    nms <- names(d)
    name_bytes <- nchar(nms, type = "bytes")
    if (any(name_bytes > 1000L) || anyDuplicated(nms)) {
      if (!skip_invalid) {
        stop("File ", basename(f), " has invalid column names (oversized or duplicated).",
             call. = FALSE)
      }
      warning("Skipping ", basename(f),
              ": invalid column names (oversized or duplicated).", call. = FALSE)
      skipped <- c(skipped, basename(f))
      next
    }
    out[[f]] <- d
  }
  if (length(out) == 0L) {
    stop("No valid files could be read from ", dir, call. = FALSE)
  }
  if (length(skipped)) {
    message("Skipped invalid file(s): ", paste(skipped, collapse = ", "))
  }
  names(out) <- vapply(out, function(d) d[["sample_id"]][1], character(1))
  if (combine) {
    out <- dplyr::bind_rows(out)
  }
  out
}
