#' Per-FOV quality control summary
#'
#' Computes per-FOV QC metrics from annotated cell data (mirroring
#' `mIF_Tzone_05_QC_Visual.R`): total cell number, DC and T cell counts and
#' their ratios, the mean cell count per sample and a low-cell-number flag
#' (a FOV is flagged `"LowCellNum"` when it contains fewer than
#' `mean_cells / low_cell_factor` cells).
#'
#' @param data A data frame with FOV, cell-type and sample columns (output of
#'   [assign_fov()] applied to annotated data).
#' @param fov_col Name of the FOV column.
#' @param celltype_col Name of the cell-type column.
#' @param sample_col Name of the sample column.
#' @param dc_types Character vector of cell types counted as DC. If `NULL`,
#'   types whose name contains `"DC"` are used (e.g. `"cDC"`,
#'   `"cDC(Mature/Active)"`).
#' @param tcell_types Character vector of cell types counted as T cells. If
#'   `NULL`, types starting with `"CD4"` or `"CD8"` are used.
#' @param low_cell_factor Denominator of the low-cell-number threshold
#'   (default 4, i.e. a FOV with fewer than one quarter of the sample mean is
#'   flagged).
#' @param max_dc_ratio Maximum accepted fraction of DC cells in a FOV. FOVs
#'   above this value are flagged as `"HighDCRatio"`. The default is `0.70`;
#'   use `NULL` to disable this rule.
#'
#' @return A data frame with one row per FOV and columns `Sample`, FOV,
#'   `TotalCells`, `DC_count`, `Tcell_count`, `DC_ratio`, `Tcell_ratio`,
#'   `MeanCells`, `LowCellFlag`, `HighDCRatioFlag`, and `InitialQC`.
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
#' head(qc)
qc_fov_summary <- function(data,
                           fov_col = "FOV",
                           celltype_col = "CellType",
                           sample_col = "sample_id",
                           dc_types = NULL,
                           tcell_types = NULL,
                           low_cell_factor = 4,
                           max_dc_ratio = 0.70) {
  fc <- .find_col(data, fov_col)
  cc <- .find_col(data, celltype_col)
  sc <- .find_col(data, sample_col)
  if (is.null(fc) || is.null(cc) || is.null(sc)) {
    stop("Could not find the FOV, cell-type and/or sample columns. Pass ",
         "`fov_col`, `celltype_col` and `sample_col` explicitly.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    stop("`data` must contain at least one cell.", call. = FALSE)
  }
  if (!is.numeric(low_cell_factor) || length(low_cell_factor) != 1L ||
      !is.finite(low_cell_factor) || low_cell_factor <= 0) {
    stop("`low_cell_factor` must be one positive finite number.",
         call. = FALSE)
  }
  if (!is.null(max_dc_ratio) &&
      (!is.numeric(max_dc_ratio) || length(max_dc_ratio) != 1L ||
       !is.finite(max_dc_ratio) || max_dc_ratio < 0 || max_dc_ratio > 1)) {
    stop("`max_dc_ratio` must be NULL or one finite number in [0, 1].",
         call. = FALSE)
  }
  types <- unique(data[[cc]])
  if (is.null(dc_types)) dc_types <- grep("DC", types, value = TRUE)
  if (is.null(tcell_types)) {
    tcell_types <- types[grepl("^CD4|^CD8", types)]
  }
  if (length(dc_types) == 0L && length(tcell_types) == 0L) {
    warning("No DC or T-cell types detected; counts will be 0. Pass ",
            "`dc_types` / `tcell_types` explicitly.", call. = FALSE)
  }

  summ <- data |>
    dplyr::mutate(
      .__sample = .data[[sc]],
      .__dc = .data[[cc]] %in% dc_types,
      .__tcell = .data[[cc]] %in% tcell_types
    ) |>
    dplyr::group_by(.__sample, .data[[fc]]) |>
    dplyr::summarise(
      Sample = dplyr::first(.__sample),
      TotalCells = dplyr::n(),
      DC_count = sum(.__dc),
      Tcell_count = sum(.__tcell),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      DC_ratio = DC_count / TotalCells,
      Tcell_ratio = Tcell_count / TotalCells
    ) |>
    dplyr::group_by(Sample) |>
    dplyr::mutate(
      MeanCells = mean(TotalCells),
      LowCellFlag = ifelse(
        TotalCells < MeanCells / low_cell_factor, "LowCellNum", "OK"
      ),
      HighDCRatioFlag = if (is.null(max_dc_ratio)) {
        "NotEvaluated"
      } else {
        ifelse(.data$DC_ratio > max_dc_ratio, "HighDCRatio", "OK")
      },
      InitialQC = ifelse(
        .data$LowCellFlag == "OK" &
          .data$HighDCRatioFlag %in% c("OK", "NotEvaluated"),
        "OK",
        ifelse(
          .data$LowCellFlag != "OK" &
            .data$HighDCRatioFlag == "HighDCRatio",
          "LowCellNum;HighDCRatio",
          ifelse(.data$LowCellFlag != "OK", "LowCellNum", "HighDCRatio")
        )
      )
    ) |>
    dplyr::ungroup()

  summ[, c("Sample", fc, "TotalCells", "DC_count", "Tcell_count",
           "DC_ratio", "Tcell_ratio", "MeanCells", "LowCellFlag",
           "HighDCRatioFlag", "InitialQC"),
       drop = FALSE]
}
