#' Bar chart of nearest-neighbour cell-type ranking
#'
#' Visualises the output of [nearest_neighbor_rank()]: for each nearest cell
#' type, the number of FOVs in which it appears among the `top_n` most
#' frequent nearest neighbours. Optionally faceted by sample. Uses the
#' [ggplot2::theme_classic()] publication style.
#'
#' @param data Output of [nearest_neighbor_rank()] (or any data frame with
#'   nearest-type, rank and count columns).
#' @param top_n Keep only types whose `Rank` is `<= top_n` (default 3).
#' @param nearest_col Name of the nearest-neighbour type column.
#' @param rank_col Name of the rank column.
#' @param count_col Name of the count column (used only for validation).
#' @param sample_col Name of the sample column (used when `by_sample = TRUE`).
#' @param by_sample Logical; if `TRUE`, bars are split/faceted by sample.
#' @param fill_colour Fill colour for the aggregated bars.
#' @param palette Optional character vector of colours (when `by_sample =
#'   TRUE`).
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
#' nb <- find_neighbors(fov, radius = 60)
#' rank <- nearest_neighbor_rank(nb)
#' plot_nearest_rank(rank, top_n = 3)
plot_nearest_rank <- function(data,
                              top_n = 3,
                              nearest_col = "nearestNGB",
                              rank_col = "Rank",
                              count_col = "Count",
                              sample_col = "sample_id",
                              by_sample = FALSE,
                              fill_colour = "#4DBBD5",
                              palette = NULL,
                              title = NULL) {
  nc <- .find_col(data, nearest_col)
  rc <- .find_col(data, rank_col)
  cc <- .find_col(data, count_col)
  if (is.null(nc) || is.null(rc) || is.null(cc)) {
    stop("Could not find the nearest-type, rank and/or count columns.",
         call. = FALSE)
  }
  d <- data[!is.na(data[[rc]]) & data[[rc]] <= top_n, , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No rows with Rank <= top_n after filtering.", call. = FALSE)
  }

  y_lab <- paste0("FOVs where type is Top", top_n)
  if (!by_sample) {
    agg <- d |>
      dplyr::count(.data[[nc]], name = "FOV_count") |>
      dplyr::arrange(dplyr::desc(FOV_count))
    p <- ggplot2::ggplot(
      agg,
      ggplot2::aes(x = stats::reorder(.data[[nc]], FOV_count), y = FOV_count)
    ) +
      ggplot2::geom_col(fill = fill_colour, width = 0.7)
  } else {
    sc <- .find_col(data, c(sample_col, "Sample", "sample_id"))
    if (is.null(sc)) {
      stop("Could not find the sample column. Pass `sample_col`.",
           call. = FALSE)
    }
    agg <- d |>
      dplyr::count(.data[[sc]], .data[[nc]], name = "FOV_count") |>
      dplyr::arrange(dplyr::desc(FOV_count))
    n_samp <- dplyr::n_distinct(agg[[sc]])
    if (is.null(palette)) palette <- .default_palette(n_samp)
    p <- ggplot2::ggplot(
      agg,
      ggplot2::aes(
        x = stats::reorder(.data[[nc]], FOV_count),
        y = FOV_count,
        fill = .data[[sc]]
      )
    ) +
      ggplot2::geom_col(position = "dodge", width = 0.7) +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::labs(fill = sc)
  }

  p <- p +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::labs(x = NULL, y = y_lab, title = title)
  p
}