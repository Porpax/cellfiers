#' Remove outliers with the IQR rule
#'
#' Drops rows whose `value_col` lies outside
#' `[Q1 - k * IQR, Q3 + k * IQR]`, computed either globally or per group
#' (e.g. per cell type). This is the outlier cleaning step used before box
#' plots in `mIF_Tzone_05_QC_Visual.R`.
#'
#' @param data A data frame.
#' @param value_col Name of the numeric column to clean.
#' @param group_col Optional name of a grouping column; bounds are computed
#'   within each group.
#' @param k IQR multiplier (default 1.5).
#'
#' @return The data frame with outlier rows removed.
#' @export
#'
#' @examples
#' df <- data.frame(
#'   group = rep(c("A", "B"), each = 10),
#'   value = c(rnorm(10, 5, 1), rnorm(10, 50, 10))
#' )
#' clean <- remove_outliers_iqr(df, value_col = "value", group_col = "group")
#' nrow(clean)
remove_outliers_iqr <- function(data,
                                value_col,
                                group_col = NULL,
                                k = 1.5) {
  vc <- .find_col(data, value_col)
  if (is.null(vc)) {
    stop("Could not find the value column '", value_col, "'.", call. = FALSE)
  }
  if (!is.numeric(data[[vc]])) {
    stop("`value_col` must be numeric.", call. = FALSE)
  }
  if (is.null(group_col)) {
    qs <- stats::quantile(data[[vc]], c(0.25, 0.75), na.rm = TRUE)
    iqr <- qs[2] - qs[1]
    lo <- qs[1] - k * iqr
    hi <- qs[2] + k * iqr
    data[!is.na(data[[vc]]) & data[[vc]] >= lo & data[[vc]] <= hi, , drop = FALSE]
  } else {
    gc <- .find_col(data, group_col)
    if (is.null(gc)) {
      stop("Could not find the group column '", group_col, "'.", call. = FALSE)
    }
    bounds <- data |>
      dplyr::group_by(.data[[gc]]) |>
      dplyr::summarise(
        .__q1 = stats::quantile(.data[[vc]], 0.25, na.rm = TRUE),
        .__q3 = stats::quantile(.data[[vc]], 0.75, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        .__lo = .__q1 - k * (.__q3 - .__q1),
        .__hi = .__q3 + k * (.__q3 - .__q1)
      )
    data |>
      dplyr::left_join(bounds, by = gc) |>
      dplyr::filter(
        !is.na(.data[[vc]]) &
          .data[[vc]] >= .__lo & .data[[vc]] <= .__hi
      ) |>
      dplyr::select(-.__q1, -.__q3, -.__lo, -.__hi)
  }
}