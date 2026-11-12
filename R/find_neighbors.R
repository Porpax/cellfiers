# encoding: UTF-8
#' Find neighbouring cells within a radius (per FOV)
#'
#' For every cell, searches for cells of the same FOV within `radius`
#' micrometres using a KD-tree ([RANN::nn2()]). This mirrors the neighbour
#' analysis of `mIF_Tzone_03_FOV_cell_Distance.R`.
#'
#' For each centre cell the result contains:
#' \itemize{
#'   \item one count column per cell type (number of neighbours of that type
#'     within the radius);
#'   \item `n_neighbors`: total number of neighbours within the radius;
#'   \item `nearestNGB`: cell type(s) of the closest neighbour(s) (ties are
#'     comma-joined), or `"None"` when the cell has no neighbours.
#' }
#'
#' @param data A data frame with FOV, cell-type and centroid columns (output
#'   of [assign_fov()] applied to annotated data).
#' @param radius Search radius in micrometres.
#'   Radii above 30 micrometres are allowed but trigger a warning because
#'   runtime and neighbour density increase and biological interpretation may
#'   become less local.
#' @param celltype_col Name of the cell-type column.
#' @param fov_col Name of the FOV column.
#' @param x,y Names of the centroid columns. Defaults to QuPath-style
#'   `"Centroid X µm"` / `"Centroid Y µm"`.
#' @param id_col Optional column with unique cell ids (e.g. `cell_id`). If
#'   `NULL`, a synthetic id is used internally.
#' @param max_k Maximum number of neighbours requested from the KD-tree
#'   (performance cap).
#' @param include_self Logical; if `TRUE`, the centre cell itself counts as a
#'   neighbour (default `FALSE`).
#'
#' @return The input data frame with neighbour-count columns, `n_neighbors`
#'   and `nearestNGB` added. Row order is preserved.
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
#' head(nb[, c("cell_id", "n_neighbors", "nearestNGB")])
find_neighbors <- function(data,
                           radius = 20,
                           celltype_col = "CellType",
                           fov_col = "FOV",
                           x = NULL,
                           y = NULL,
                           id_col = NULL,
                           max_k = 1000,
                           include_self = FALSE) {
  cc <- .find_col(data, celltype_col)
  fc <- .find_col(data, fov_col)
  xc <- .find_col(data, c(x, "Centroid X um", "Centroid X", "x", "X"))
  yc <- .find_col(data, c(y, "Centroid Y um", "Centroid Y", "y", "Y"))
  if (is.null(cc) || is.null(fc) || is.null(xc) || is.null(yc)) {
    stop("Could not find the cell-type, FOV and/or centroid columns. Pass ",
         "`celltype_col`, `fov_col`, `x` and `y` explicitly.", call. = FALSE)
  }
  if (nrow(data) == 0L) stop("`data` must contain at least one cell.", call. = FALSE)
  if (!is.numeric(radius) || length(radius) != 1L ||
      !is.finite(radius) || radius <= 0) {
    stop("`radius` must be one positive finite number.", call. = FALSE)
  }
  if (radius > 30) {
    warning(
      "`radius` is greater than 30 um. Computation and memory use may ",
      "increase, and the biological meaning of a local neighbourhood may ",
      "weaken. Confirm that this radius is appropriate for your assay.",
      call. = FALSE
    )
  }
  if (!is.numeric(max_k) || length(max_k) != 1L ||
      !is.finite(max_k) || max_k < 1) {
    stop("`max_k` must be one positive finite number.", call. = FALSE)
  }
  max_k <- as.integer(max_k)
  if (anyNA(data[[fc]]) || anyNA(data[[cc]])) {
    stop("FOV and cell-type columns must not contain missing values.",
         call. = FALSE)
  }
  if (!is.numeric(data[[xc]]) || !is.numeric(data[[yc]]) ||
      any(!is.finite(data[[xc]])) || any(!is.finite(data[[yc]]))) {
    stop("Centroid coordinates must be finite numeric values.", call. = FALSE)
  }

  types <- unique(data[[cc]])
  new_cols <- setdiff(types, names(data))
  collisions <- intersect(types, names(data))
  if (length(collisions)) {
    warning("Cell-type names collide with existing columns: ",
            paste(collisions, collapse = ", "),
            "; existing columns will be overwritten by neighbour counts.",
            call. = FALSE)
  }

  idc <- .find_col(data, c(id_col, "cell_id"))
  if (is.null(idc)) {
    data[[".__id"]] <- paste0("cell_", seq_len(nrow(data)))
  } else {
    data[[".__id"]] <- as.character(data[[idc]])
  }
  data[[".__rowid"]] <- seq_len(nrow(data))

  fovs <- unique(data[[fc]])
  out_list <- vector("list", length(fovs))
  capped <- FALSE
  for (fi in seq_along(fovs)) {
    f <- fovs[fi]
    idx <- which(data[[fc]] == f)
    d <- data[idx, , drop = FALSE]
    coords <- as.matrix(d[, c(xc, yc)])
    rownames(coords) <- d[[".__id"]]
    k <- min(nrow(coords), max_k)
    nn <- RANN::nn2(coords, query = coords, k = k, radius = radius)
    if (k < nrow(coords) && any(nn$nn.idx[, k] > 0L)) capped <- TRUE

    counts <- matrix(0L, nrow = nrow(d), ncol = length(types),
                     dimnames = list(NULL, types))
    n_total <- integer(nrow(d))
    nearest <- rep(NA_character_, nrow(d))
    nb_type <- d[[cc]]

    for (j in seq_len(nrow(d))) {
      idxx <- nn$nn.idx[j, ]
      dst <- nn$nn.dists[j, ]
      if (include_self) {
        valid <- !is.na(dst) & dst <= radius & idxx > 0
      } else {
        valid <- !is.na(dst) & dst <= radius & idxx > 0 & idxx != j
      }
      if (!any(valid)) next
      nb_idx <- idxx[valid]
      tt <- table(factor(nb_type[nb_idx], levels = types))
      counts[j, ] <- as.integer(tt)
      n_total[j] <- length(nb_idx)
      md <- min(dst[valid])
      nearest[j] <- paste(unique(nb_type[nb_idx][dst[valid] == md]),
                          collapse = ",")
    }

    for (t in types) d[[t]] <- counts[, t]
    d[["n_neighbors"]] <- n_total
    d[["nearestNGB"]] <- nearest
    d[["nearestNGB"]][is.na(d[["nearestNGB"]])] <- "None"
    out_list[[fi]] <- d
  }

  res <- dplyr::bind_rows(out_list)
  res <- res[order(res[[".__rowid"]]), , drop = FALSE]
  res[[".__rowid"]] <- NULL
  res[[".__id"]] <- NULL
  if (capped) {
    warning("At least one cell reached `max_k`; neighbour counts may be ",
            "truncated. Increase `max_k` and rerun.", call. = FALSE)
  }
  res
}
