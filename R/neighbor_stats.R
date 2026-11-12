#' Per-FOV neighbour statistics for centre-cell types
#'
#' Two summary functions built on top of [find_neighbors()]:
#'
#' - `neighbor_proportions()`: keeps centre cells of the requested types and,
#'   for each of them, computes the total number of neighbours and the
#'   proportion of neighbours of each cell type (mirroring the per-DC
#'   computation of `mIF_Tzone_04_statistic.R`).
#' - `neighbor_fov_summary()`: aggregates neighbour counts and proportions per
#'   FOV (mean count and mean proportion per neighbour cell type, plus the
#'   number of centre cells per FOV), in long (default) or wide format.
#'
#' @param data Output of [find_neighbors()].
#' @param center_types Cell types to use as centre cells (e.g.
#'   `c("cDC", "cDC(Mature/Active)")`). If `NULL`, all cell types are used.
#' @param celltype_col Name of the cell-type column.
#' @param fov_col Name of the FOV column.
#' @param sample_col Name of the sample column.
#' @param wide For `neighbor_fov_summary()`, return wide format (one row per
#'   FOV) instead of long format.
#'
#' @return
#' - `neighbor_proportions()`: the centre-cell data frame with
#'   `total_neighbors` and `Pro_<type>` columns.
#' - `neighbor_fov_summary()`: a data frame with columns `Sample`, FOV,
#'   `CellType`, `mean_count`, `mean_prop`, `n_center` (long), or one row per
#'   FOV with `mean_count_<type>` / `mean_prop_<type>` columns (wide).
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
#' fov <- assign_fov(annotated)
#' nb <- find_neighbors(fov, radius = 60)
#' dc <- neighbor_proportions(nb, center_types = c("cDC", "cDC(Mature/Active)"))
#' summ <- neighbor_fov_summary(nb, center_types = c("cDC", "cDC(Mature/Active)"))
#' head(summ)
neighbor_proportions <- function(data,
                                 center_types = NULL,
                                 celltype_col = "CellType") {
  cc <- .find_col(data, celltype_col)
  if (is.null(cc)) {
    stop("Could not find the cell-type column.", call. = FALSE)
  }
  types <- unique(data[[cc]])
  neighbor_cols <- setdiff(intersect(types, names(data)), cc)
  if (length(neighbor_cols) == 0L) {
    stop("No neighbour-count columns found. Run find_neighbors() first.",
         call. = FALSE)
  }
  if (is.null(center_types)) center_types <- types
  d <- data[data[[cc]] %in% center_types, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No centre cells of the requested types.", call. = FALSE)
  }
  total <- rowSums(d[, neighbor_cols, drop = FALSE])
  d[["total_neighbors"]] <- total
  d <- d[total > 0, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No centre cells have neighbours within the search radius.",
         call. = FALSE)
  }
  for (t in neighbor_cols) {
    d[[paste0("Pro_", t)]] <- d[[t]] / d[["total_neighbors"]]
  }
  d
}

#' @rdname neighbor_proportions
#' @export
neighbor_fov_summary <- function(data,
                                 center_types = NULL,
                                 celltype_col = "CellType",
                                 fov_col = "FOV",
                                 sample_col = "sample_id",
                                 wide = FALSE) {
  cc <- .find_col(data, celltype_col)
  fc <- .find_col(data, fov_col)
  sc <- .find_col(data, sample_col)
  if (is.null(cc) || is.null(fc) || is.null(sc)) {
    stop("Could not find the cell-type, FOV and/or sample columns.",
         call. = FALSE)
  }
  types <- unique(data[[cc]])
  neighbor_cols <- setdiff(intersect(types, names(data)), cc)
  if (length(neighbor_cols) == 0L) {
    stop("No neighbour-count columns found. Run find_neighbors() first.",
         call. = FALSE)
  }
  if (is.null(center_types)) center_types <- types
  d <- data[data[[cc]] %in% center_types, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No centre cells of the requested types.", call. = FALSE)
  }
  total <- rowSums(d[, neighbor_cols, drop = FALSE])
  d[["total_neighbors"]] <- total
  d <- d[total > 0, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No centre cells have neighbours within the search radius.",
         call. = FALSE)
  }

  # the centre-cell type column is dropped; neighbour types become CellType
  d[[cc]] <- NULL
  long <- d |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(neighbor_cols),
      names_to = "CellType",
      values_to = "count"
    ) |>
    dplyr::mutate(prop = count / total_neighbors) |>
    dplyr::group_by(.data[[fc]], CellType) |>
    dplyr::summarise(
      Sample = dplyr::first(.data[[sc]]),
      mean_count = mean(count),
      mean_prop = mean(prop),
      n_center = dplyr::n(),
      .groups = "drop"
    )
  long <- long[, c("Sample", fc, "CellType", "mean_count", "mean_prop",
                   "n_center"), drop = FALSE]

  if (!wide) {
    return(long)
  }
  wide_out <- long |>
    tidyr::pivot_wider(
      id_cols = c("Sample", dplyr::all_of(fc), "n_center"),
      names_from = CellType,
      values_from = c(mean_count, mean_prop),
      names_sep = "_"
    )
  wide_out
}

#' Summarise FOV neighbour statistics at the biological-sample level
#'
#' Averages FOV-level neighbour statistics within each sample and neighbour
#' cell type. This creates one observation per biological sample and cell type,
#' which is generally more appropriate than treating FOVs as independent
#' replicates in between-group tests.
#'
#' @param data Long-format output of [neighbor_fov_summary()].
#' @param sample_col Name of the sample column.
#' @param celltype_col Name of the neighbour cell-type column.
#'
#' @return A data frame with `Sample`, `CellType`, `mean_count`, `mean_prop`,
#'   `n_fov`, and `n_center` columns.
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
#' nb <- analyze_neighbors(assign_fov(annotated), radius = 60)
#' neighbor_sample_summary(nb$fov_summary)
neighbor_sample_summary <- function(data,
                                    sample_col = "Sample",
                                    celltype_col = "CellType") {
  sc <- .find_col(data, sample_col)
  cc <- .find_col(data, celltype_col)
  required <- c("mean_count", "mean_prop", "n_center")
  if (is.null(sc) || is.null(cc) || !all(required %in% names(data))) {
    stop("`data` must be a long-format FOV neighbour summary.",
         call. = FALSE)
  }
  if (nrow(data) == 0L) stop("`data` must contain at least one row.", call. = FALSE)

  out <- data |>
    dplyr::mutate(
      .__sample = .data[[sc]],
      .__celltype = .data[[cc]]
    ) |>
    dplyr::group_by(.data$.__sample, .data$.__celltype) |>
    dplyr::summarise(
      mean_count = mean(mean_count, na.rm = TRUE),
      mean_prop = mean(mean_prop, na.rm = TRUE),
      n_fov = dplyr::n(),
      n_center = sum(n_center),
      .groups = "drop"
    )
  names(out)[names(out) == ".__sample"] <- "Sample"
  names(out)[names(out) == ".__celltype"] <- "CellType"
  out
}
