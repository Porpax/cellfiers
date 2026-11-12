# encoding: UTF-8
#' Plot the FOV spatial distribution
#'
#' Mirrors the FOV distribution figure of the original workflow
#' (`mIF_Tzone_02_GridArea_Process.R`): cell centroids coloured by FOV with
#' `geom_point(size = 0.4, alpha = 0.6)`, a turbo colour ramp from the viridis
#' package, `theme_classic()` without grid lines, fixed aspect ratio, and axis
#' labels "Centroid X (µm)" / "Centroid Y (µm)". The title is
#' "FOV distribution - <sample>".
#'
#' @param data A data frame with FOV labels and centroid coordinates (output
#'   of [assign_fov()]).
#' @param fov_col Name of the FOV column.
#' @param sample_col Name of the sample column (used for faceting/filtering).
#' @param sample_id Optional sample id; if given, only that sample is plotted.
#' @param x,y Names of the centroid columns. Defaults to QuPath-style
#'   `"Centroid X µm"` / `"Centroid Y µm"`.
#' @param point_size Point size (default 0.4, as in the original).
#' @param alpha Point transparency (default 0.6, as in the original).
#' @param palette Optional character vector of colours (one per FOV). If
#'   `NULL`, a turbo palette is generated with `viridis::viridis(30,
#'   option = "turbo")` (recycled and shuffled when > 30 FOVs).
#' @param show_legend Logical; if `TRUE` show the FOV legend (default `FALSE`
#'   as in the original).
#' @param title Optional plot title (defaults to "FOV distribution -
#'   <sample>").
#' @param x_lab,y_lab Axis labels (defaults match the original figure).
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
#' plot_fov_distribution(fov)
#' plot_fov_distribution(fov, sample_id = "C1001")
plot_fov <- function(data,
                                  fov_col = "FOV",
                                  sample_col = "sample_id",
                                  sample_id = NULL,
                                  x = NULL,
                                  y = NULL,
                                  point_size = 0.4,
                                  alpha = 0.6,
                                  palette = NULL,
                                  show_legend = FALSE,
                                  title = NULL,
                                  x_lab = "Centroid X (\u00b5m)",
                                  y_lab = "Centroid Y (\u00b5m)") {
  fc <- .find_col(data, fov_col)
  xc <- .find_col(data, c(x, "Centroid X um", "Centroid X", "x", "X"))
  yc <- .find_col(data, c(y, "Centroid Y um", "Centroid Y", "y", "Y"))
  sc <- .find_col(data, sample_col)
  if (is.null(fc) || is.null(xc) || is.null(yc)) {
    stop("Could not find the FOV and/or centroid columns. Pass `fov_col`, ",
         "`x` and `y` explicitly.", call. = FALSE)
  }

  if (!is.null(sample_id)) {
    if (is.null(sc)) {
      stop("`sample_id` was given but no sample column was found.",
           call. = FALSE)
    }
    data <- data[data[[sc]] == sample_id, , drop = FALSE]
    if (nrow(data) == 0L) {
      stop("No cells found for sample_id '", sample_id, "'.", call. = FALSE)
    }
  }

  fovs <- unique(data[[fc]])
  n_fov <- length(fovs)
  if (is.null(palette)) {
    base <- .turbo_palette()
    if (n_fov <= 30L) {
      palette <- base[seq_len(n_fov)]
    } else {
      palette <- rep(base, ceiling(n_fov / 30L))[seq_len(n_fov)]
      set.seed(123)
      palette <- sample(palette)
    }
  }

  if (is.null(title) && !is.null(sample_id)) {
    title <- paste("FOV distribution -", sample_id)
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(
    x = .data[[xc]], y = .data[[yc]], color = .data[[fc]]
  )) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = title, x = x_lab, y = y_lab) +
    ggplot2::coord_fixed() +
    ggplot2::scale_color_manual(values = palette)

  if (!is.null(sc) && is.null(sample_id) && dplyr::n_distinct(data[[sc]]) > 1) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", sc)))
  }
  p
}

#' @rdname plot_fov
#' @export
plot_fov_distribution <- plot_fov
