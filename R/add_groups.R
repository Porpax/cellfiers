#' Add experimental group information to a data frame
#'
#' Left-joins a sample-to-group mapping (such as
#' `Docu/SLN_mIF_Group_pCRpNR.txt`, columns `Sample` and `Group`) onto cell or
#' FOV level data.
#'
#' @param data A data frame with a sample column.
#' @param groups A data frame with sample ids and group labels.
#' @param sample_col Name of the sample column in `data`.
#' @param sample_id_col Name of the sample-id column in `groups` (defaults to
#'   `"Sample"`).
#' @param group_col Name of the output group column (defaults to `"Group"`).
#'
#' @return The input data frame with an extra `group_col` column.
#' @export
#'
#' @examples
#' data(mif_groups)
#' data(mif_cells)
#' with_group <- add_groups(mif_cells, mif_groups)
#' head(with_group[, c("sample_id", "Group")])
add_groups <- function(data,
                       groups,
                       sample_col = "sample_id",
                       sample_id_col = "Sample",
                       group_col = "Group") {
  sc <- .find_col(data, sample_col)
  gs <- .find_col(groups, c(sample_id_col, "Sample"))
  gg <- .find_col(groups, c(group_col, "Group"))
  if (is.null(sc)) {
    stop("Could not find the sample column in `data`.", call. = FALSE)
  }
  if (is.null(gs) || is.null(gg)) {
    stop("Could not find the sample-id and/or group column in `groups`.",
         call. = FALSE)
  }
  if (group_col %in% names(data)) {
    warning("Column '", group_col, "' already exists and will be overwritten.",
            call. = FALSE)
  }
  groups_use <- groups[, c(gs, gg), drop = FALSE]
  names(groups_use) <- c(".__sample_key", ".__group")
  res <- data |>
    dplyr::mutate(.__sample_key = .data[[sc]]) |>
    dplyr::left_join(groups_use, by = ".__sample_key")
  res[[group_col]] <- res[[".__group"]]
  res[, setdiff(names(res), c(".__sample_key", ".__group")), drop = FALSE]
}