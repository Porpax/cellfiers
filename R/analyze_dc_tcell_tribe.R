#' Analyse DC-centred target-T neighbourhoods (compatibility wrapper)
#'
#' This legacy name is retained for existing scripts. New analyses should use
#' [target_neighbors()], which accepts arbitrary centre and target
#' cell types.
#'
#' @inheritParams target_neighbors
#' @param center_types Centre cell types.
#' @param target_t_types Target cell types.
#' @param other_t_types Optional other cell types used for labelling.
#'
#' @return A list containing the legacy `dc_level` element plus `fov_level`,
#'   `edges`, `settings`, and `weighted_columns`.
#' @examples
#' cells <- data.frame(
#'   sample_id = "S1", FOV = "F1",
#'   CellType = c("cDC", "CD8Tn(Primed)"), x = c(0, 10), y = c(0, 0)
#' )
#' analyze_dc_tcell_tribe(cells)$dc_level
#' @export
analyze_dc_tcell_tribe <- function(
    data,
    radius = 20,
    center_types = c("cDC", "cDC(Mature/Active)"),
    target_t_types = c(
      "CD8Tn(Primed)", "CD8Teff(Active)",
      "CD4Tn(Primed)", "CD4Teff(Active)"
    ),
    other_t_types = c(
      "CD8Tn(Quiescent)", "CD8Teff(Resting)",
      "CD4Tn(Quiescent)", "CD4Teff(Resting)"
    ),
    shared_weight = c("inverse_degree", "none"),
    sample_col = "sample_id",
    celltype_col = "CellType",
    fov_col = "FOV",
    x = NULL,
    y = NULL,
    grid_size = 200,
    max_k = 1000L) {
  result <- target_neighbors(
    data = data,
    center_types = center_types,
    target_types = target_t_types,
    radius = radius,
    other_types = other_t_types,
    shared_weight = match.arg(shared_weight),
    sample_col = sample_col,
    celltype_col = celltype_col,
    fov_col = fov_col,
    x = x, y = y, grid_size = grid_size, max_k = max_k
  )
  result$dc_level <- result$center_level
  result
}
