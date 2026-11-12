#' Plot the original-workflow DC-count density
#'
#' Reproduces the density figures in `mIF_Tzone_05_QC_Visual.R` without
#' changing their geometry, colours, typography, or labels.
#'
#' @param data A data frame containing a numeric DC-count column.
#' @param metric Name of the DC-count column.
#' @param filtered Logical; use the original title for the central 90 percent
#'   filtered data when `TRUE`.
#'
#' @return A [ggplot2::ggplot()] object. Save at 8 by 4 inches to reproduce
#'   the original output dimensions.
#' @export
#'
#' @examples
#' qc <- data.frame(DC_count = c(0, 1, 2, 2, 3, 5, 8))
#' plot_dc_count_distribution(qc)
plot_dc_counts <- function(data,
                                       metric = "DC_count",
                                       filtered = FALSE) {
  if (!metric %in% names(data) || !is.numeric(data[[metric]])) {
    stop("`metric` must name a numeric column in `data`.", call. = FALSE)
  }
  if (nrow(data) < 2L || length(unique(stats::na.omit(data[[metric]]))) < 2L) {
    stop("At least two distinct non-missing values are required.",
         call. = FALSE)
  }
  plot_title <- if (isTRUE(filtered)) {
    "Distribution of DC Counts 90% in all FOVs"
  } else {
    "Distribution of DC Counts in all FOVs"
  }

  ggplot2::ggplot(data, ggplot2::aes(x = .data[[metric]])) +
    ggplot2::geom_density(fill = "steelblue", alpha = 0.5, color = "black") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "DC Count", y = "Density", title = plot_title) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black", size = 10),
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

#' @rdname plot_dc_counts
#' @export
plot_dc_count_distribution <- plot_dc_counts
