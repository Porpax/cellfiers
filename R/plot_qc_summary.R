#' Plot QC metrics by sample
#'
#' Box plot (with optional jitter points) of one or more per-FOV QC metrics,
#' faceted by metric, in the original workflow style: `theme_classic()`,
#' black axis text (size 10), bold axis titles (size 11).
#'
#' @param qc The QC summary returned by [qc_fov_summary()].
#' @param metrics Character vector of numeric QC metrics to plot.
#' @param sample_col Name of the sample column in `qc`.
#' @param color_by Optional column (e.g. `LowCellFlag`) used to colour the
#'   points.
#' @param jitter Logical; if `TRUE`, individual FOV values are shown as
#'   jittered points.
#' @param palette Optional character vector of colours.
#' @param title Optional plot title.
#'
#' @return A [ggplot2::ggplot()] object.
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
#' plot_qc_summary(qc)
plot_qc_summary <- function(qc,
                            metrics = c("TotalCells", "DC_count",
                                        "Tcell_count", "DC_ratio",
                                        "Tcell_ratio"),
                            sample_col = "Sample",
                            color_by = NULL,
                            jitter = TRUE,
                            palette = NULL,
                            title = NULL) {
  sc <- .find_col(qc, sample_col)
  if (is.null(sc)) {
    stop("Could not find the sample column in `qc`.", call. = FALSE)
  }
  metrics <- intersect(metrics, names(qc))
  if (length(metrics) == 0L) {
    stop("None of the requested `metrics` exist in `qc`.", call. = FALSE)
  }
  cb <- .find_col(qc, c(color_by, "LowCellFlag"))

  d <- qc |>
    dplyr::select(dplyr::all_of(c(sc, metrics, cb))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(metrics),
      names_to = "metric",
      values_to = "value"
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data[[sc]], y = value, fill = .data[[sc]]
  )) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(x = NULL, y = "Value", title = title)

  if (jitter) {
    if (is.null(cb)) {
      p <- p + ggplot2::geom_jitter(
        width = 0.2, size = 0.5, alpha = 0.4
      )
    } else {
      p <- p + ggplot2::geom_jitter(
        ggplot2::aes(color = .data[[cb]]),
        width = 0.2, size = 0.5
      )
    }
  }

  n_samp <- dplyr::n_distinct(d[[sc]])
  if (is.null(palette)) palette <- .default_palette(n_samp)
  p <- p + ggplot2::scale_fill_manual(values = palette)
  if (jitter && !is.null(cb)) {
    n_flag <- dplyr::n_distinct(d[[cb]])
    p <- p + ggplot2::scale_color_manual(values = .default_palette(n_flag))
  }
  p
}

#' Plot the density of a QC metric (original workflow style)
#'
#' Density plot of a per-FOV QC metric (e.g. `DC_count`), mirroring the QC
#' density figures of `mIF_Tzone_05_QC_Visual.R`: steelblue fill,
#' `alpha = 0.5`, black outline, `theme_classic()` with black axis text and
#' bold axis titles. By default all samples are pooled into one density;
#' pass `by_sample = TRUE` to colour by sample instead.
#'
#' @param qc The QC summary returned by [qc_fov_summary()].
#' @param metric Numeric QC metric to plot.
#' @param sample_col Name of the sample column in `qc`.
#' @param by_sample Logical; if `TRUE`, densities are coloured by sample.
#' @param cutoffs Optional numeric vector; vertical dashed lines drawn at
#'   these values (e.g. the quantile boundaries).
#' @param palette Optional character vector of colours (when
#'   `by_sample = TRUE`).
#' @param title Optional plot title.
#'
#' @return A [ggplot2::ggplot()] object.
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
#' plot_metric_density(qc, metric = "DC_count")
plot_metric_density <- function(qc,
                                metric = "DC_count",
                                sample_col = "Sample",
                                by_sample = FALSE,
                                cutoffs = NULL,
                                palette = NULL,
                                title = NULL) {
  if (!by_sample && is.null(cutoffs) && is.null(title) &&
      identical(metric, "DC_count")) {
    return(plot_dc_count_distribution(qc, metric = metric, filtered = FALSE))
  }
  if (!metric %in% names(qc)) {
    stop("Metric '", metric, "' not found in `qc`.", call. = FALSE)
  }
  sc <- .find_col(qc, sample_col)

  if (by_sample && !is.null(sc)) {
    p <- ggplot2::ggplot(qc, ggplot2::aes(
      x = .data[[metric]], fill = .data[[sc]]
    )) +
      ggplot2::geom_density(alpha = 0.5, color = "black")
    n_samp <- dplyr::n_distinct(qc[[sc]])
    if (is.null(palette)) palette <- .default_palette(n_samp)
    p <- p + ggplot2::scale_fill_manual(values = palette)
  } else {
    p <- ggplot2::ggplot(qc, ggplot2::aes(x = .data[[metric]])) +
      ggplot2::geom_density(fill = "steelblue", alpha = 0.5, color = "black")
  }

  if (!is.null(cutoffs)) {
    cut_df <- data.frame(xintercept = cutoffs)
    p <- p + ggplot2::geom_vline(
      data = cut_df, ggplot2::aes(xintercept = xintercept),
      linetype = "dashed", color = "grey30"
    )
  }

  p <- p +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    ) +
    ggplot2::labs(x = metric, y = "Density", title = title)
  p
}
