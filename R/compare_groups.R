#' Compare groups with a statistical test
#'
#' Runs a Wilcoxon rank-sum (default) or Student t test between groups for a
#' numeric variable, either globally or separately for each level of a
#' `by` column (e.g. cell type). All pairwise group combinations are tested.
#' Observations are assumed to be independent. For imaging studies with
#' multiple FOVs per biological sample, first aggregate with
#' [neighbor_sample_summary()] to avoid treating FOVs as biological replicates.
#'
#' @param data A data frame with a numeric value column, a group column and
#'   optionally a `by` column.
#' @param value_col Name of the numeric column to compare.
#' @param group_col Name of the group column.
#' @param by Optional column whose levels are compared separately (e.g.
#'   `"CellType"`).
#' @param method Test to use: `"wilcox"` (Wilcoxon rank-sum / Mann-Whitney)
#'   or `"t"` (Student t test).
#' @param paired Logical; if `TRUE`, a paired test is used.
#'
#' @return A data frame with one row per comparison: `group1`, `group2`,
#'   `statistic`, `p.value`, `n1`, `n2`, `method`, plus the `by` column when
#'   provided.
#' @export
#'
#' @examples
#' df <- data.frame(
#'   Group = rep(c("R", "NR"), each = 10),
#'   value = c(rnorm(10, 10, 2), rnorm(10, 14, 2))
#' )
#' compare_groups(df, value_col = "value", group_col = "Group")
compare_groups <- function(data,
                           value_col,
                           group_col,
                           by = NULL,
                           method = c("wilcox", "t"),
                           paired = FALSE) {
  vc <- .find_col(data, value_col)
  gc <- .find_col(data, group_col)
  method <- match.arg(method)
  if (is.null(vc)) {
    stop("Could not find the value column '", value_col, "'.", call. = FALSE)
  }
  if (is.null(gc)) {
    stop("Could not find the group column '", group_col, "'.", call. = FALSE)
  }
  if (!is.numeric(data[[vc]])) {
    stop("`value_col` must be numeric.", call. = FALSE)
  }
  groups <- unique(data[[gc]])
  groups <- groups[!is.na(groups)]
  if (length(groups) < 2L) {
    stop("Need at least two groups in `group_col`.", call. = FALSE)
  }

  .run_test <- function(x, y) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    if (length(x) < 2L || length(y) < 2L) {
      return(list(p.value = NA_real_, statistic = NA_real_,
                  n1 = length(x), n2 = length(y)))
    }
    if (method == "wilcox") {
      res <- suppressWarnings(
        stats::wilcox.test(x, y, paired = paired, exact = FALSE)
      )
    } else {
      res <- stats::t.test(x, y, paired = paired)
    }
    list(p.value = res$p.value, statistic = unname(res$statistic),
         n1 = length(x), n2 = length(y))
  }

  .compare_one <- function(sub) {
    pairs <- utils::combn(groups, 2)
    rows <- lapply(seq_len(ncol(pairs)), function(i) {
      g1 <- pairs[1, i]
      g2 <- pairs[2, i]
      t <- .run_test(sub[[vc]][sub[[gc]] == g1], sub[[vc]][sub[[gc]] == g2])
      data.frame(
        group1 = g1, group2 = g2,
        statistic = t$statistic, p.value = t$p.value,
        n1 = t$n1, n2 = t$n2, method = method,
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(rows)
  }

  if (is.null(by)) {
    return(.compare_one(data))
  }
  bc <- .find_col(data, by)
  if (is.null(bc)) {
    stop("Could not find the `by` column '", by, "'.", call. = FALSE)
  }
  levels_by <- unique(data[[bc]])
  rows <- lapply(levels_by, function(lv) {
    sub <- data[data[[bc]] == lv, , drop = FALSE]
    r <- .compare_one(sub)
    if (nrow(r)) r[[bc]] <- lv
    r
  })
  dplyr::bind_rows(rows)
}
