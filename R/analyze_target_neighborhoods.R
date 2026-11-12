# encoding: UTF-8
#' Analyse user-defined centre--target neighbourhoods
#'
#' Searches around any user-selected centre cell type and quantifies any set
#' of target cell types. A target shared by `m` centre cells is weighted by
#' `1/m` by default, preventing the same target cell from being counted in
#' full for every nearby centre.
#'
#' @param data One biological sample of annotated, FOV-assigned cell data.
#' @param center_types Character vector of cell types used as centres.
#' @param target_types Character vector of target cell types used to calculate
#'   richness, raw abundance and weighted abundance.
#' @param radius Search radius in micrometres. Values above 30 are allowed but
#'   trigger a warning about computation and biological interpretation.
#' @param other_types Optional non-target types used to label an
#'   `Other-selected` neighbourhood when no requested target is present.
#' @param shared_weight Either `"inverse_degree"` (default, `1/m`) or
#'   `"none"` (every centre--target edge has weight 1).
#' @param sample_col,celltype_col,fov_col Column names.
#' @param x,y Optional centroid column names.
#' @param grid_size FOV side length in micrometres, used for map coordinates.
#' @param max_k Maximum KD-tree neighbours requested per centre.
#'
#' @return A list with `center_level`, `fov_level`, `edges`, `settings`, and
#'   the target-composition column names in `weighted_columns`.
#' @examples
#' cells <- data.frame(
#'   sample_id = "S1", FOV = "F1",
#'   CellType = c("DC", "DC", "CD4 T", "CD8 T"),
#'   x = c(0, 10, 5, 15), y = c(0, 0, 0, 0)
#' )
#' result <- target_neighbors(cells, "DC", c("CD4 T", "CD8 T"), radius = 20)
#' result$fov_level
#' @export
target_neighbors <- function(
    data,
    center_types,
    target_types,
    radius = 20,
    other_types = NULL,
    shared_weight = c("inverse_degree", "none"),
    sample_col = "sample_id",
    celltype_col = "CellType",
    fov_col = "FOV",
    x = NULL,
    y = NULL,
    grid_size = 200,
    max_k = 1000L) {
  shared_weight <- match.arg(shared_weight)
  sc <- .find_col(data, sample_col)
  cc <- .find_col(data, celltype_col)
  fc <- .find_col(data, fov_col)
  xc <- .find_col(data, c(x, "Centroid X um", "Centroid X", "x", "X"))
  yc <- .find_col(data, c(y, "Centroid Y um", "Centroid Y", "y", "Y"))
  if (any(vapply(list(sc, cc, fc, xc, yc), is.null, logical(1)))) {
    stop("Could not find sample, cell type, FOV, or centroid columns.",
         call. = FALSE)
  }
  sample_ids <- unique(as.character(data[[sc]]))
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  if (length(sample_ids) != 1L) {
    stop("`data` must contain exactly one non-missing sample.", call. = FALSE)
  }
  if (!length(center_types) || !length(target_types) ||
      anyNA(center_types) || anyNA(target_types)) {
    stop("`center_types` and `target_types` must be non-empty.",
         call. = FALSE)
  }
  if (length(intersect(center_types, target_types))) {
    warning("Some cell types occur in both `center_types` and `target_types`. ",
            "This is allowed, but verify that it matches your question.",
            call. = FALSE)
  }
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
  if (!is.numeric(grid_size) || length(grid_size) != 1L ||
      !is.finite(grid_size) || grid_size <= 0) {
    stop("`grid_size` must be one positive finite number.", call. = FALSE)
  }
  if (!is.numeric(data[[xc]]) || !is.numeric(data[[yc]]) ||
      any(!is.finite(data[[xc]])) || any(!is.finite(data[[yc]]))) {
    stop("Centroid coordinates must be finite numeric values.", call. = FALSE)
  }

  center_rows <- which(as.character(data[[cc]]) %in% center_types)
  if (!length(center_rows)) {
    stop("No requested centre cells were found in sample ", sample_ids,
         call. = FALSE)
  }
  coordinates <- as.matrix(data[, c(xc, yc)])
  k_use <- min(as.integer(max_k), nrow(data))
  nearest <- RANN::nn2(
    data = coordinates,
    query = coordinates[center_rows, , drop = FALSE],
    k = k_use, radius = radius, searchtype = "standard"
  )
  if (k_use < nrow(data) && any(nearest$nn.dists[, k_use] <= radius)) {
    warning("At least one centre reached `max_k`; increase it and rerun.",
            call. = FALSE)
  }
  neighbor_index <- lapply(seq_along(center_rows), function(i) {
    idx <- nearest$nn.idx[i, ]
    distance <- nearest$nn.dists[i, ]
    valid <- idx > 0L & is.finite(distance) & distance <= radius &
      idx != center_rows[i]
    idx[valid]
  })

  edges <- dplyr::bind_rows(lapply(seq_along(center_rows), function(i) {
    idx <- neighbor_index[[i]]
    keep <- as.character(data[[cc]][idx]) %in% target_types
    if (!any(keep)) return(NULL)
    data.frame(
      Sample = sample_ids,
      CenterIndex = i,
      CenterCellRow = center_rows[i],
      CenterType = as.character(data[[cc]][center_rows[i]]),
      FOV = as.character(data[[fc]][center_rows[i]]),
      TargetCellRow = idx[keep],
      TargetType = as.character(data[[cc]][idx[keep]]),
      stringsAsFactors = FALSE
    )
  }))
  if (nrow(edges)) {
    edges <- edges |>
      dplyr::group_by(.data$TargetCellRow) |>
      dplyr::mutate(SharedByM = dplyr::n_distinct(.data$CenterIndex)) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        EdgeWeight = if (shared_weight == "inverse_degree") {
          1 / .data$SharedByM
        } else {
          1
        }
      )
  } else {
    edges <- data.frame(
      Sample = character(), CenterIndex = integer(),
      CenterCellRow = integer(), CenterType = character(), FOV = character(),
      TargetCellRow = integer(), TargetType = character(),
      SharedByM = integer(), EdgeWeight = numeric()
    )
  }

  center_level <- dplyr::bind_rows(lapply(seq_along(center_rows), function(i) {
    idx <- neighbor_index[[i]]
    types <- as.character(data[[cc]][idx])
    raw_counts <- vapply(
      target_types, function(z) sum(types == z, na.rm = TRUE), integer(1)
    )
    center_edges <- edges[edges$CenterIndex == i, , drop = FALSE]
    weighted_counts <- vapply(target_types, function(z) {
      sum(center_edges$EdgeWeight[center_edges$TargetType == z], na.rm = TRUE)
    }, numeric(1))
    richness <- sum(raw_counts > 0L)
    other_count <- if (length(other_types)) {
      sum(types %in% other_types)
    } else {
      0L
    }
    centre_neighbours <- sum(types %in% center_types)
    neighborhood_class <- if (richness > 0L) {
      paste0("Target-L", richness)
    } else if (other_count > 0L) {
      "Other-selected"
    } else if (centre_neighbours > 0L) {
      "Centre-only"
    } else if (!length(idx)) {
      "Isolated"
    } else {
      "Other-cell"
    }
    out <- data.frame(
      Sample = sample_ids,
      CenterIndex = i,
      CenterCellRow = center_rows[i],
      CenterType = as.character(data[[cc]][center_rows[i]]),
      FOV = as.character(data[[fc]][center_rows[i]]),
      Radius_um = radius,
      TotalNeighbors = length(idx),
      TargetRichness = richness,
      RawTargetAbundance = sum(raw_counts),
      WeightedTargetAbundance = sum(weighted_counts),
      SharedTargetCount = sum(center_edges$SharedByM > 1L),
      OtherSelectedCount = other_count,
      CentreNeighborCount = centre_neighbours,
      NeighborhoodClass = neighborhood_class,
      stringsAsFactors = FALSE
    )
    for (z in seq_along(target_types)) {
      out[[paste0("Raw_", target_types[z])]] <- raw_counts[z]
      out[[paste0("Weighted_", target_types[z])]] <- weighted_counts[z]
    }
    out
  }))

  weighted_columns <- paste0("Weighted_", target_types)
  all_fovs <- data.frame(
    Sample = sample_ids,
    FOV = as.character(data[[fc]]),
    FOV_center_X = (floor(data[[xc]] / grid_size) + 0.5) * grid_size,
    FOV_center_Y = (floor(data[[yc]] / grid_size) + 0.5) * grid_size,
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$Sample, .data$FOV) |>
    dplyr::summarise(
      FOV_center_X = dplyr::first(.data$FOV_center_X),
      FOV_center_Y = dplyr::first(.data$FOV_center_Y),
      .groups = "drop"
    )

  fov_rows <- lapply(split(seq_len(nrow(center_level)), center_level$FOV),
                     function(ii) {
    z <- center_level[ii, , drop = FALSE]
    out <- data.frame(
      Sample = z$Sample[1], FOV = z$FOV[1], N_Center = nrow(z),
      MeanTargetRichness = mean(z$TargetRichness),
      MeanRawTargetAbundancePerCenter = mean(z$RawTargetAbundance),
      MeanWeightedTargetAbundancePerCenter = mean(z$WeightedTargetAbundance),
      FractionWithAnyTarget = mean(z$TargetRichness >= 1),
      FractionHighRichness = mean(z$TargetRichness >=
                                    max(1, length(target_types) - 1)),
      N_TargetNeighborhood = sum(z$TargetRichness >= 1),
      N_OtherSelectedNeighborhood =
        sum(z$NeighborhoodClass == "Other-selected"),
      N_CentreOnlyNeighborhood = sum(z$NeighborhoodClass == "Centre-only"),
      stringsAsFactors = FALSE
    )
    for (col in weighted_columns) out[[col]] <- mean(z[[col]])
    out
  })
  centre_fovs <- dplyr::bind_rows(fov_rows)
  fov_level <- all_fovs |>
    dplyr::left_join(centre_fovs, by = c("Sample", "FOV"))
  numeric_zero <- c(
    "N_Center", "MeanTargetRichness", "MeanRawTargetAbundancePerCenter",
    "MeanWeightedTargetAbundancePerCenter", "FractionWithAnyTarget",
    "FractionHighRichness", "N_TargetNeighborhood",
    "N_OtherSelectedNeighborhood", "N_CentreOnlyNeighborhood",
    weighted_columns
  )
  for (col in numeric_zero) fov_level[[col]][is.na(fov_level[[col]])] <- 0
  fov_level$FOVClass <- ifelse(
    fov_level$N_Center == 0, "No-centre",
    ifelse(
      fov_level$N_TargetNeighborhood > 0, "Target",
      ifelse(
        fov_level$N_OtherSelectedNeighborhood > 0, "Other-selected",
        ifelse(fov_level$N_CentreOnlyNeighborhood > 0,
               "Centre-only", "Other/isolated")
      )
    )
  )

  list(
    center_level = center_level,
    fov_level = fov_level,
    edges = edges,
    weighted_columns = weighted_columns,
    settings = list(
      sample = sample_ids, center_types = center_types,
      target_types = target_types, other_types = other_types,
      radius = radius, shared_weight = shared_weight,
      grid_size = grid_size, max_k = max_k
    )
  )
}

#' @rdname target_neighbors
#' @export
analyze_target_neighborhoods <- target_neighbors
