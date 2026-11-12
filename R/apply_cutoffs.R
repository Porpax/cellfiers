#' Apply fluorescence channel cutoffs to call markers positive/negative
#'
#' For each marker in `marker_map`, cells with a channel intensity greater
#' than or equal to the sample-specific cutoff are labelled positive, others
#' negative. Cutoffs are stored in a data frame such as
#' `Docu/SLN_mIF_6channel_cutoff.txt` with one row per sample and one column
#' per marker.
#'
#' @param data A data frame with QuPath cell data (must contain a sample
#'   column and the intensity columns named by `marker_map`).
#' @param cutoffs A data frame with sample ids and one numeric cutoff column
#'   per marker. The first column (or `cutoff_id_col`) holds the sample ids.
#' @param marker_map A named character vector mapping marker names to
#'   intensity column names, e.g.
#'   `c(CD8 = "Cell: Opal 780 mean", CD3 = "Cell: Opal 480 mean")`.
#' @param sample_col Name of the sample column in `data`.
#' @param cutoff_id_col Name of the sample-id column in `cutoffs`. Defaults
#'   to the first column of `cutoffs`.
#' @param pos Character used for positive calls (default `"+"`).
#' @param neg Character used for negative calls (default `"-"`).
#'
#' @return The input data frame with one extra column per marker holding
#'   `pos`/`neg` (or `NA` when no cutoff is available for the sample, or when
#'   the intensity is missing).
#' @export
#'
#' @examples
#' data(mif_cells)
#' data(mif_cutoffs)
#' marker_map <- c(
#'   CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
#'   CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
#'   PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
#' )
#' cells_pos <- apply_cutoffs(mif_cells, mif_cutoffs, marker_map)
apply_cutoffs <- function(data,
                          cutoffs,
                          marker_map,
                          sample_col = "sample_id",
                          cutoff_id_col = NULL,
                          pos = "+",
                          neg = "-") {
  sc <- .find_col(data, sample_col)
  if (is.null(sc)) {
    stop("Could not find the sample column. Pass `sample_col` explicitly.",
         call. = FALSE)
  }
  if (is.null(cutoff_id_col)) cutoff_id_col <- names(cutoffs)[1]
  if (!cutoff_id_col %in% names(cutoffs)) {
    stop("cutoff_id_col '", cutoff_id_col, "' not found in `cutoffs`.",
         call. = FALSE)
  }
  if (is.null(names(marker_map)) || any(!nzchar(names(marker_map)))) {
    stop("`marker_map` must be a named character vector.", call. = FALSE)
  }

  samples <- unique(data[[sc]])
  for (marker in names(marker_map)) {
    int_col <- marker_map[[marker]]
    if (!int_col %in% names(data)) {
      stop("Intensity column '", int_col, "' (for marker ", marker,
           ") not found in `data`.", call. = FALSE)
    }
    if (!marker %in% names(cutoffs)) {
      stop("Cutoff column '", marker, "' not found in `cutoffs`.",
           call. = FALSE)
    }
    intensity <- suppressWarnings(as.numeric(data[[int_col]]))
    if (any(is.na(intensity) & !is.na(data[[int_col]]))) {
      warning("Non-numeric values found in '", int_col,
              "'; they are treated as NA.", call. = FALSE)
    }
    cutoff_vec <- rep(NA_real_, nrow(data))
    for (s in samples) {
      hit <- which(cutoffs[[cutoff_id_col]] == s)
      if (length(hit) == 0L) {
        warning("No cutoff row for sample '", s, "'; marker ", marker,
                " set to NA.", call. = FALSE)
        next
      }
      cutoff_vec[data[[sc]] == s] <- as.numeric(cutoffs[[marker]][hit[1]])
    }
    data[[marker]] <- ifelse(intensity >= cutoff_vec, pos, neg)
  }
  data
}