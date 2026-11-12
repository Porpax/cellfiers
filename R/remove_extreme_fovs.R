#' Remove FOVs with extreme values of a QC metric
#'
#' Drops FOVs whose value of `metric` (e.g. `DC_count`) lies outside the
#' central `probs` quantile interval, mirroring the "remove extreme DC
#' counts" step of `mIF_Tzone_05_QC_Visual.R`. When cell-level `data` is
#' supplied, both the filtered QC summary and the filtered cell data are
#' returned.
#'
#' @param qc The QC summary returned by [qc_fov_summary()].
#' @param metric Name of the numeric QC metric column to filter on.
#' @param probs Quantile probabilities defining the kept interval
#'   (default `c(0.05, 0.95)`).
#' @param data Optional cell-level data with an FOV column; if given, it is
#'   filtered to the retained FOVs.
#' @param fov_col Name of the FOV column.
#'
#' @return If `data` is `NULL`, the filtered QC summary. Otherwise a list with
#'   elements `qc` (filtered summary) and `data` (filtered cell data).
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
#' qc <- qc_fov_summary(fov)
#' qc_f <- remove_extreme_fovs(qc, metric = "DC_count")
remove_extreme_fovs <- function(qc,
                                metric = "DC_count",
                                probs = c(0.05, 0.95),
                                data = NULL,
                                fov_col = "FOV") {
  fc <- .find_col(qc, fov_col)
  if (is.null(fc)) {
    stop("Could not find the FOV column in `qc`.", call. = FALSE)
  }
  if (!metric %in% names(qc)) {
    stop("Metric '", metric, "' not found in `qc`.", call. = FALSE)
  }
  if (length(probs) != 2L || probs[1] < 0 || probs[2] > 1 || probs[1] >= probs[2]) {
    stop("`probs` must be two increasing probabilities in [0, 1].",
         call. = FALSE)
  }
  q <- stats::quantile(qc[[metric]], probs = probs, na.rm = TRUE)
  keep_fovs <- qc[[fc]][!is.na(qc[[metric]]) &
                           qc[[metric]] >= q[1] & qc[[metric]] <= q[2]]
  out_qc <- qc[qc[[fc]] %in% keep_fovs, , drop = FALSE]

  if (is.null(data)) {
    return(out_qc)
  }
  d_fc <- .find_col(data, fov_col)
  if (is.null(d_fc)) {
    stop("Could not find the FOV column in `data`.", call. = FALSE)
  }
  out_data <- data[data[[d_fc]] %in% keep_fovs, , drop = FALSE]
  list(qc = out_qc, data = out_data)
}