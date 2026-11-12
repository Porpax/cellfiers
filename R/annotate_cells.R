#' Assign cell types from a classification rule table
#'
#' Each cell is matched against the rows of a rule table (such as
#' `Docu/SLN_mIF_anno_rule.txt`). A rule row lists, for each marker, a
#' required status (`"+"`, `"-"` or `NA` = "don't care"). The first rule
#' whose non-NA markers all agree with the cell's marker status wins. Cells
#' matching no rule are labelled `unmatched` (default `"Other"`).
#'
#' Marker status columns in `data` may hold `"+"`/`"-"`, `1`/`0` or
#' `TRUE`/`FALSE`; they are normalised internally.
#'
#' @param data A data frame with marker status columns (e.g. produced by
#'   [apply_cutoffs()]).
#' @param rules A data frame whose first column (or `rule_col`) holds the
#'   cell type labels and whose remaining columns hold marker statuses.
#' @param markers Optional character vector of marker columns to use. Defaults
#'   to the columns of `rules` other than `rule_col` that also exist in
#'   `data`.
#' @param celltype_col Name of the output column holding the assigned cell
#'   type. Defaults to `"CellType"`.
#' @param rule_col Name of the cell-type column in `rules`. Defaults to the
#'   first column of `rules`.
#' @param unmatched Label for cells that match no rule (default `"Other"`).
#'
#' @return The input data frame with an extra `celltype_col` column.
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
#' pos <- apply_cutoffs(mif_cells, mif_cutoffs, marker_map)
#' annotated <- annotate_cells(pos, mif_rules)
annotate_cells <- function(data,
                           rules,
                           markers = NULL,
                           celltype_col = "CellType",
                           rule_col = NULL,
                           unmatched = "Other") {
  if (is.null(rule_col)) rule_col <- names(rules)[1]
  if (!rule_col %in% names(rules)) {
    stop("Cell-type column '", rule_col, "' not found in `rules`.", call. = FALSE)
  }
  if (is.null(markers)) {
    markers <- setdiff(intersect(names(rules), names(data)), rule_col)
  }
  if (length(markers) == 0L) {
    stop("No marker columns shared between `rules` and `data`.", call. = FALSE)
  }
  miss <- setdiff(markers, names(data))
  if (length(miss)) {
    stop("Marker column(s) not found in `data`: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  n <- nrow(data)
  nrule <- nrow(rules)
  cell_mat <- vapply(markers, function(m) .normalize_marker(data[[m]]),
                     character(n))
  if (is.null(dim(cell_mat))) cell_mat <- matrix(cell_mat, nrow = n)
  rule_mat <- vapply(markers, function(m) .normalize_marker(rules[[m]]),
                     character(nrule))
  if (is.null(dim(rule_mat))) rule_mat <- matrix(rule_mat, nrow = nrule)

  out <- rep(unmatched, n)
  unmatched_idx <- seq_len(n)
  for (r in seq_len(nrule)) {
    if (!length(unmatched_idx)) break
    ok <- rep(TRUE, length(unmatched_idx))
    for (k in seq_along(markers)) {
      rv <- rule_mat[r, k]
      if (!is.na(rv)) {
        cv <- cell_mat[unmatched_idx, k]
        ok <- ok & (!is.na(cv)) & (cv == rv)
      }
    }
    if (any(ok)) {
      hit <- unmatched_idx[ok]
      out[hit] <- rules[[rule_col]][r]
      unmatched_idx <- unmatched_idx[!ok]
    }
  }
  data[[celltype_col]] <- out
  data
}