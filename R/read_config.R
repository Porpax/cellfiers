#' Read configuration tables used by the annotation pipeline
#'
#' Convenience wrappers around [readr::read_tsv()] for the tab-delimited
#' configuration files used by the mIF workflow:
#'
#' - `read_cutoffs()`: per-sample fluorescence channel cutoffs (first column
#'   sample id, remaining columns marker names).
#' - `read_rules()`: cell-type classification rules (first column cell type,
#'   remaining columns marker statuses `"+"`/`"-"`/`NA`).
#' - `read_groups()`: sample-to-group mapping (columns `Sample`, `Group`).
#'
#' @param path Path to a tab-delimited file.
#' @param ... Further arguments passed to [readr::read_tsv()].
#'
#' @return A [tibble::tibble()].
#' @export
#'
#' @examples
#' \dontrun{
#' cutoffs <- read_cutoffs("Docu/SLN_mIF_6channel_cutoff.txt")
#' rules <- read_rules("Docu/SLN_mIF_anno_rule.txt")
#' groups <- read_groups("Docu/SLN_mIF_Group_pCRpNR.txt")
#' }
read_cutoffs <- function(path, ...) {
  readr::read_tsv(path, show_col_types = FALSE, ...)
}

#' @rdname read_cutoffs
#' @export
read_rules <- function(path, ...) {
  readr::read_tsv(path, show_col_types = FALSE, ...)
}

#' @rdname read_cutoffs
#' @export
read_groups <- function(path, ...) {
  readr::read_tsv(path, show_col_types = FALSE, ...)
}