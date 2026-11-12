#' Rasterise cells and run field-of-view quality control
#'
#' Runs the Phase 2 workflow in its intended order: assign FOVs, flag FOVs
#' below the per-sample mean cell count divided by `low_cell_factor`, flag
#' FOVs above `max_dc_ratio`, and then remove FOVs outside the requested
#' quantiles of a QC metric. Quantiles are estimated only from FOVs that
#' passed the low-cell and DC-ratio rules.
#'
#' @inheritParams assign_fov
#' @inheritParams qc_fov_summary
#' @param metric Numeric QC metric used for the extreme-FOV filter.
#' @param probs Two probabilities defining the retained quantile interval.
#'
#' @return A named list containing `data_fov`, `qc`, `qc_pass`, `qc_final`,
#'   `data_final`, and `metric_quantiles`.
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
#' result <- run_fov_qc(annotated, probs = c(0, 1))
#' names(result)
run_fov_qc <- function(data,
                       sample_col = "sample_id",
                       x = NULL,
                       y = NULL,
                       grid_size = 200,
                       min_distance = 500,
                       fov_col = "FOV",
                       celltype_col = "CellType",
                       dc_types = NULL,
                       tcell_types = NULL,
                       low_cell_factor = 4,
                       max_dc_ratio = 0.70,
                       metric = "DC_count",
                       probs = c(0.05, 0.95)) {
  data_fov <- assign_fov(
    data = data,
    sample_col = sample_col,
    x = x,
    y = y,
    grid_size = grid_size,
    min_distance = min_distance,
    fov_col = fov_col
  )
  qc <- qc_fov_summary(
    data = data_fov,
    fov_col = fov_col,
    celltype_col = celltype_col,
    sample_col = sample_col,
    dc_types = dc_types,
    tcell_types = tcell_types,
    low_cell_factor = low_cell_factor,
    max_dc_ratio = max_dc_ratio
  )
  qc_pass <- qc[qc$InitialQC == "OK", , drop = FALSE]
  if (nrow(qc_pass) == 0L) {
    stop("No FOV passed the low-cell-number filter.", call. = FALSE)
  }
  if (!metric %in% names(qc_pass) || !is.numeric(qc_pass[[metric]])) {
    stop("`metric` must name a numeric column in the QC summary.",
         call. = FALSE)
  }
  data_pass <- data_fov[data_fov[[fov_col]] %in% qc_pass[[fov_col]], , drop = FALSE]
  metric_quantiles <- stats::quantile(
    qc_pass[[metric]], probs = probs, na.rm = TRUE, names = TRUE
  )
  final <- remove_extreme_fovs(
    qc_pass, metric = metric, probs = probs, data = data_pass,
    fov_col = fov_col
  )
  qc$ExtremeMetricFlag <- ifelse(
    qc$InitialQC != "OK",
    "NotEvaluated",
    ifelse(
      !is.na(qc[[metric]]) & qc[[metric]] >= metric_quantiles[[1]] &
        qc[[metric]] <= metric_quantiles[[2]],
      "OK", "ExtremeMetric"
    )
  )
  qc$FinalQC <- ifelse(
    qc$InitialQC == "OK" & qc$ExtremeMetricFlag == "OK", "OK", "Exclude"
  )
  list(
    data_fov = data_fov,
    qc = qc,
    qc_pass = qc_pass,
    qc_final = final$qc,
    data_final = final$data,
    metric_quantiles = metric_quantiles
  )
}
