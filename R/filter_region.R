#' Filter cells by parent region
#'
#' Keeps (or drops) cells whose parent annotation matches user-selected
#' regions. Region names are not restricted to `Tregion` or `SLN_Tregion`;
#' any QuPath annotation name can be selected by exact, fixed-substring or
#' regular-expression matching.
#'
#' @param data A data frame with QuPath cell data.
#' @param parent_col Name of the parent/annotation column. If NULL, the
#'   function looks for a column named "Parent" (or a normalised match).
#' @param regions Character vector of parent values to keep (or drop).
#' @param keep Logical; if TRUE keep the matching rows, if FALSE drop them.
#' @param match Matching mode: `"exact"`, `"fixed"` substring, or `"regex"`.
#' @param ignore_case Logical; ignore case for `fixed` or `regex` matching.
#'
#' @return A filtered data frame.
#' @export
#'
#' @examples
#' data(mif_cells)
#' tzone <- filter_region(mif_cells, regions = "SLN_Tregion")
#' tumour <- filter_region(mif_cells, regions = "Tumou?r|Tumor",
#'                         match = "regex", ignore_case = TRUE)
filter_region <- function(data,
                          parent_col = NULL,
                          regions,
                          keep = TRUE,
                          match = c("exact", "fixed", "regex"),
                          ignore_case = FALSE) {
  match <- match.arg(match)
  pc <- .find_col(data, c(parent_col, "Parent"))
  if (is.null(pc)) {
    stop("Could not find the parent column. Pass parent_col explicitly.",
         call. = FALSE)
  }
  if (!length(regions) || anyNA(regions)) {
    stop("`regions` must contain at least one non-missing region pattern.",
         call. = FALSE)
  }
  vals <- as.character(data[[pc]])
  if (match == "exact") {
    if (isTRUE(ignore_case)) {
      flag <- tolower(vals) %in% tolower(regions)
    } else {
      flag <- vals %in% regions
    }
  } else {
    matcher <- if (match == "fixed") {
      function(pattern) grepl(pattern, vals, fixed = TRUE,
                              ignore.case = ignore_case)
    } else {
      function(pattern) grepl(pattern, vals, perl = TRUE,
                              ignore.case = ignore_case)
    }
    flag <- Reduce(`|`, lapply(regions, matcher))
  }
  flag[is.na(flag)] <- FALSE
  if (!keep) flag <- !flag
  data[flag, , drop = FALSE]
}
