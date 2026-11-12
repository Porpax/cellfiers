#' Radar chart of cell-type composition by group (fmsb style)
#'
#' Draws a group-comparison radar chart with [fmsb::radarchart()], following
#' the reference fmsb code of the original workflow: axes start at the top and
#' are connected clockwise, polygons are filled with a semi-transparent colour
#' (`alpha = 0.1`), borders are solid (`plwd = 1`), grid lines are grey
#' (`cglwd = 0.8`), axis labels are black, and the legend is drawn on the
#' right with `pch = 20`. Cell types listed in `exclude` (by default
#' `"Other"`) are removed before computing proportions, and (for
#' `value = "proportion"`) values are shown as percentages with `seg = 4`
#' grid segments.
#'
#' Because the chart uses base graphics, it is drawn either to `file` (when
#' given) or to the current graphics device.
#'
#' @param data A data frame with cell-type, sample and group columns (e.g.
#'   annotated cells after [add_groups()]).
#' @param celltype_col Name of the cell-type column.
#' @param sample_col Name of the sample column.
#' @param group_col Name of the group column.
#' @param value What to average per group: `"proportion"` (default, cell-type
#'   proportion within each sample, in percent) or `"count"` (mean cell
#'   count).
#' @param exclude Cell types to remove before analysis (default `"Other"`).
#'   Use `NULL` to keep every cell type.
#' @param percent Logical; if `TRUE` (default) and `value = "proportion"`,
#'   values are multiplied by 100 and the axis labels show percentages.
#' @param levels Optional ordering of cell types around the radar (defaults to
#'   the order of appearance).
#' @param palette Optional character vector of colours, one per group. Use a
#'   named vector to fix colours to groups, e.g.
#'   `c(R = "#E64B35", NR = "#4DBBD5")`. Defaults to the reference palette
#'   `c("#1B9E77", "#7570B3", "#D95F02")`.
#' @param file Optional output file (`.png`, `.jpg`, `.pdf`, `.bmp`). If
#'   `NULL`, the chart is drawn on the current graphics device.
#' @param width,height Figure size in centimetres (for raster output).
#' @param res Resolution in dpi for raster output (default 600, as in the
#'   reference).
#' @param family Font family (default `"serif"`, as in the reference).
#' @param seg Number of grid segments (default 4).
#' @param plwd Polygon border line width (default 1).
#' @param cglwd Grid line width (default 0.8).
#' @param vlcex Vertex (cell-type) label size (default 0.8).
#' @param legend_cex Legend text size (default 1.2).
#' @param title Optional plot title.
#' @param ... Further arguments passed to [fmsb::radarchart()].
#'
#' @return Invisibly, the data frame used by `radarchart` (rows: `max`, `min`,
#'   then one row per group; columns: cell types).
#' @export
#'
#' @examples
#' data(mif_cells)
#' data(mif_cutoffs)
#' data(mif_rules)
#' data(mif_groups)
#' marker_map <- c(
#'   CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
#'   CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
#'   PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
#' )
#' annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
#' annotated <- add_groups(annotated, mif_groups)
#' \dontrun{
#' plot_celltype_radar(annotated)
#' plot_celltype_radar(annotated, file = "radar.png")
#' }
plot_celltype_radar <- function(data,
                                celltype_col = "CellType",
                                sample_col = "sample_id",
                                group_col = "Group",
                                value = c("proportion", "count"),
                                exclude = "Other",
                                percent = TRUE,
                                levels = NULL,
                                palette = NULL,
                                file = NULL,
                                width = 15,
                                height = 15,
                                res = 600,
                                family = "serif",
                                seg = 4,
                                plwd = 1,
                                cglwd = 0.8,
                                vlcex = 0.8,
                                legend_cex = 1.2,
                                title = NULL,
                                ...) {
  cc <- .find_col(data, celltype_col)
  sc <- .find_col(data, sample_col)
  gc <- .find_col(data, group_col)
  value <- match.arg(value)
  if (is.null(cc) || is.null(sc) || is.null(gc)) {
    stop("Could not find the cell-type, sample and/or group columns.",
         call. = FALSE)
  }
  if (!is.null(exclude)) {
    data <- data[!data[[cc]] %in% exclude, , drop = FALSE]
  }
  if (nrow(data) == 0L) {
    stop("No cells left after excluding the requested cell types.",
         call. = FALSE)
  }
  if (is.null(levels)) levels <- unique(data[[cc]])

  d <- data |>
    dplyr::mutate(
      .__sample = .data[[sc]],
      .__group = .data[[gc]],
      .__cc = factor(.data[[cc]], levels = levels)
    ) |>
    dplyr::count(.__sample, .__group, .__cc, .drop = FALSE) |>
    dplyr::group_by(.__sample) |>
    dplyr::mutate(.__total = sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(value = if (value == "proportion") n / .__total else n)

  agg <- d |>
    dplyr::group_by(.__group, .__cc) |>
    dplyr::summarise(value = mean(value), .groups = "drop") |>
    dplyr::mutate(.__cc = factor(.__cc, levels = levels))

  if (value == "proportion" && percent) {
    agg <- agg |>
      dplyr::mutate(value = value * 100)
  }

  groups <- unique(as.character(agg$.__group))
  n_groups <- length(groups)
  if (is.null(palette)) {
    ref_pal <- c("#1B9E77", "#7570B3", "#D95F02")
    if (n_groups <= length(ref_pal)) {
      palette <- ref_pal[seq_len(n_groups)]
    } else {
      palette <- .default_palette(n_groups)
    }
  }
  if (is.null(names(palette))) {
    palette <- stats::setNames(palette, groups)
  }

  # wide matrix: one row per group, one column per cell type
  wide <- agg |>
    dplyr::select(.__group, .__cc, value) |>
    tidyr::pivot_wider(names_from = .__cc, values_from = value) |>
    as.data.frame()
  rownames(wide) <- as.character(wide$.__group)
  wide$.__group <- NULL
  wide <- wide[, as.character(levels), drop = FALSE]

  maxv <- max(agg$value, na.rm = TRUE)
  step <- .nice_step(maxv / seg)
  chart_max <- step * seg
  chart_min <- 0
  df <- rbind(
    rep(chart_max, ncol(wide)),
    rep(chart_min, ncol(wide)),
    as.matrix(wide)
  )
  rownames(df) <- c("max", "min", groups)
  colnames(df) <- colnames(wide)
  df <- as.data.frame(df)

  caxislabels <- seq(chart_min, chart_max, by = step)

  if (!is.null(file)) {
    ext <- tolower(sub(".*[.]", "", file))
    switch(ext,
      png = grDevices::png(file, width = width, height = height,
                           units = "cm", res = res, family = family),
      jpg = ,
      jpeg = grDevices::jpeg(file, width = width, height = height,
                             units = "cm", res = res, family = family),
      bmp = grDevices::bmp(file, width = width, height = height,
                           units = "cm", res = res, family = family),
      pdf = grDevices::pdf(file, width = width / 2.54, height = height / 2.54,
                           family = family),
      stop("Unsupported file extension: ", ext, call. = FALSE)
    )
    on.exit(grDevices::dev.off())
  }

  graphics::par(mar = c(0, 0, 0, 8), xpd = TRUE)
  fmsb::radarchart(
    df,
    axistype = 1,
    pcol = palette,
    pfcol = grDevices::adjustcolor(palette, alpha.f = 0.1),
    plwd = plwd,
    plty = 1,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = caxislabels,
    seg = seg,
    cglwd = cglwd,
    vlcex = vlcex,
    ...
  )
  if (!is.null(title)) {
    graphics::title(main = title)
  }
  graphics::legend(
    x = 1.2, y = 1,
    legend = groups,
    bty = "n", pch = 20, col = palette,
    text.col = "black", cex = legend_cex, pt.cex = 3
  )
  invisible(df)
}