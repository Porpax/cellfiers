# encoding: UTF-8
#' Assign fields of view (FOVs) by rasterising tissue coordinates
#'
#' Reproduces the FOV rasterisation step of the original workflow
#' (`mIF_Tzone_02_GridArea_Process.R`). Cells are placed on a square grid of
#' `grid_size` micrometres. Grid cells are grouped into spatially separated
#' tissue regions by hierarchical clustering of the grid centres (single
#' linkage, cut at `min_distance`). Within each tissue region, grids are
#' re-numbered from the minimum x/y grid origin (row-major) and each cell is labelled
#' `"<sample>_TR<tregion>_FOV<fov>"`.
#'
#' @param data A data frame with sample, x and y columns (e.g. annotated cell
#'   data from [run_annotation()]).
#' @param sample_col Name of the sample column.
#' @param x,y Names of the centroid columns. Defaults to QuPath-style
#'   `"Centroid X µm"` / `"Centroid Y µm"` (or normalised matches).
#' @param grid_size Grid side length in micrometres.
#' @param min_distance Minimum distance (µm) between two grids for them to be
#'   considered separate tissue regions.
#' @param fov_col Name of the output FOV column.
#'
#' @return The input data frame with an extra `fov_col` column. Row order is
#'   preserved.
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
#' annotated_fov <- assign_fov(annotated)
#' head(annotated_fov$FOV)
assign_fov <- function(data,
                       sample_col = "sample_id",
                       x = NULL,
                       y = NULL,
                       grid_size = 200,
                       min_distance = 500,
                       fov_col = "FOV") {
  sc <- .find_col(data, sample_col)
  xc <- .find_col(data, c(x, "Centroid X um", "Centroid X", "x", "X"))
  yc <- .find_col(data, c(y, "Centroid Y um", "Centroid Y", "y", "Y"))
  if (is.null(sc)) {
    stop("Could not find the sample column. Pass `sample_col` explicitly.",
         call. = FALSE)
  }
  if (is.null(xc) || is.null(yc)) {
    stop("Could not find x/y centroid columns. Pass `x` and `y` explicitly.",
         call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`data` must contain at least one cell.", call. = FALSE)
  }
  if (anyNA(data[[sc]]) || any(!nzchar(as.character(data[[sc]])))) {
    stop("The sample column contains missing or empty values.", call. = FALSE)
  }
  if (!is.numeric(data[[xc]]) || !is.numeric(data[[yc]]) ||
      any(!is.finite(data[[xc]])) || any(!is.finite(data[[yc]]))) {
    stop("Centroid coordinates must be finite numeric values.", call. = FALSE)
  }
  if (grid_size <= 0 || min_distance <= 0) {
    stop("`grid_size` and `min_distance` must be positive.", call. = FALSE)
  }
  if (fov_col %in% names(data)) {
    warning("Column '", fov_col, "' already exists and will be overwritten.",
            call. = FALSE)
  }

  data[[".__rowid"]] <- seq_len(nrow(data))
  parts <- split(seq_len(nrow(data)), data[[sc]])
  out <- lapply(names(parts), function(s) {
    d <- data[parts[[s]], , drop = FALSE]
    .assign_fov_one(d, s, xc, yc, grid_size, min_distance, fov_col)
  })
  res <- dplyr::bind_rows(out)
  res <- res[order(res[[".__rowid"]]), , drop = FALSE]
  res[[".__rowid"]] <- NULL
  res
}

#' @keywords internal
.assign_fov_one <- function(d, sid, xc, yc, grid_size, min_distance, fov_col) {
  d[[".__gx"]] <- floor(d[[xc]] / grid_size)
  d[[".__gy"]] <- floor(d[[yc]] / grid_size)

  grid_summary <- d |>
    dplyr::group_by(.__gx, .__gy) |>
    dplyr::summarise(
      center_x = mean(.data[[xc]], na.rm = TRUE),
      center_y = mean(.data[[yc]], na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(grid_summary) <= 2L) {
    d[[".__tr"]] <- 1L
  } else {
    dm <- stats::dist(grid_summary[, c("center_x", "center_y")])
    hc <- stats::hclust(dm, method = "single")
    grid_summary[["cluster"]] <- stats::cutree(hc, h = min_distance)
    d <- d |>
      dplyr::left_join(
        grid_summary[, c(".__gx", ".__gy", "cluster")],
        by = c(".__gx", ".__gy")
      )
    d[[".__tr"]] <- d[["cluster"]]
    d[["cluster"]] <- NULL
  }

  d <- d |>
    dplyr::group_by(.__tr) |>
    dplyr::mutate(
      .__gx_n = .__gx - min(.__gx),
      .__gy_n = .__gy - min(.__gy),
      .__max_x = max(.__gx_n) + 1L,
      .__fov_num = .__gy_n * .__max_x + .__gx_n + 1L,
      .__fov = paste0(sid, "_TR", .__tr, "_FOV", .__fov_num)
    ) |>
    dplyr::ungroup()

  d[[fov_col]] <- d[[".__fov"]]
  drop_cols <- c(".__gx", ".__gy", ".__tr", ".__gx_n", ".__gy_n",
                 ".__max_x", ".__fov_num", ".__fov")
  d[, setdiff(names(d), drop_cols), drop = FALSE]
}
