#' cellfiers: Multiplex Immunofluorescence (mIF) Spatial Cell Analysis
#'
#' A toolbox for analysing multiplex immunofluorescence (mIF) single-cell data
#' exported from QuPath. Cells are annotated with cell types using
#' user-defined fluorescence-channel cutoffs and classification rules, then
#' summarised and visualised. Additional functions support data quality control,
#' spatial neighbour analysis, group comparisons and publication-ready figures.
#'
#' Cell annotation and display include reading QuPath
#' classification exports, filtering cells in a region of interest, calling markers positive/negative from
#' intensity cutoffs, assigning cell types from a rule table, summarising cell
#' type counts and plotting cell maps / proportions / marker distributions.
#'
#' Data quality control and cleaning include FOV
#' rasterisation, per-FOV QC metrics, low-cell-number flags, removal of
#' extreme FOVs and IQR outlier cleaning, with QC visualisations.
#'
#' Comparative analysis and visualisation include
#' KD-tree neighbour analysis, per-FOV neighbour statistics and
#' nearest-neighbour ranking, group comparisons (Wilcoxon / t test) and
#' publication-style violin/box plots with p-value brackets.
#'
#' @importFrom rlang .data
#' @keywords internal
"_PACKAGE"
