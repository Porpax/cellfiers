#' Plot marker intensity distributions (with optional cutoffs)
#'
#' Density plots of the fluorescence intensity of each marker, one panel per
#' marker, coloured by sample. When a cutoff table is supplied, a vertical
#' dashed line is drawn at each sample's cutoff. This helps to choose and
#' validate the channel thresholds used by [apply_cutoffs()].
#'
#' @param data A data frame with the intensity columns named by `marker_map`.
#' @param marker_map Named character vector mapping marker names to intensity
#'   columns (see [apply_cutoffs()]).
#' @param sample_col Name of the sample column.
#' @param cutoffs Optional cutoff data frame; if given, per-sample cutoff
#'   lines are drawn.
#' @param cutoff_id_col Name of the sample-id column in `cutoffs` (defaults to
#'   the first column).
#' @param facet_sample Logical; if `TRUE` (default) also facet by sample.
#' @param alpha Density fill transparency.
#' @param palette Optional character vector of colours.
#'
#' @return A [ggplot2::ggplot()] object.
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
#' plot_marker_distributions(mif_cells, marker_map, cutoffs = mif_cutoffs)
plot_markers <- function(data,
                                      marker_map,
                                      sample_col = "sample_id",
                                      cutoffs = NULL,
                                      cutoff_id_col = NULL,
                                      facet_sample = FALSE,
                                      alpha = 0.4,
                                      palette = NULL) {
  sc <- .find_col(data, sample_col)
  if (is.null(sc)) {
    stop("Could not find the sample column.", call. = FALSE)
  }
  if (is.null(names(marker_map)) || any(!nzchar(names(marker_map)))) {
    stop("`marker_map` must be a named character vector.", call. = FALSE)
  }
  miss <- setdiff(unname(marker_map), names(data))
  if (length(miss)) {
    stop("Intensity column(s) not found in `data`: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }

  d <- data |>
    dplyr::select(dplyr::all_of(c(sc, unname(marker_map)))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(unname(marker_map)),
      names_to = "marker_intensity",
      values_to = "intensity"
    ) |>
    dplyr::mutate(
      marker = names(marker_map)[match(marker_intensity, unname(marker_map))]
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = intensity, fill = .data[[sc]]
  )) +
    ggplot2::geom_density(alpha = alpha) +
    ggplot2::facet_wrap(~marker, scales = "free") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Mean fluorescence intensity",
      y = "Density",
      fill = sc
    )

  if (!is.null(cutoffs)) {
    if (is.null(cutoff_id_col)) cutoff_id_col <- names(cutoffs)[1]
    if (!cutoff_id_col %in% names(cutoffs)) {
      stop("cutoff_id_col '", cutoff_id_col, "' not found in `cutoffs`.",
           call. = FALSE)
    }
    lines_df <- cutoffs |>
      dplyr::select(dplyr::all_of(c(cutoff_id_col, names(marker_map)))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(names(marker_map)),
        names_to = "marker",
        values_to = "cutoff"
      ) |>
      dplyr::rename(sample_id = dplyr::all_of(cutoff_id_col))
    p <- p +
      ggplot2::geom_vline(
        data = lines_df,
        ggplot2::aes(xintercept = cutoff,
                     color = sample_id,
                     linetype = "cutoff"),
        show.legend = FALSE
      )
  }

  if (facet_sample) {
    p <- p + ggplot2::facet_grid(
      stats::as.formula(paste(sc, "~ marker")), scales = "free"
    )
  }

  n_samp <- dplyr::n_distinct(d[[sc]])
  if (is.null(palette)) palette <- .default_palette(n_samp)
  p <- p + ggplot2::scale_fill_manual(values = palette)
  p
}

#' @rdname plot_markers
#' @export
plot_marker_distributions <- plot_markers
