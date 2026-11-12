test_that("plot functions return ggplot objects", {
  skip_if_not_installed("ggplot2")
  data("mif_cells")
  data("mif_cutoffs")
  data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)

  p1 <- plot_cell_map(annotated)
  expect_s3_class(p1, "ggplot")

  p2 <- plot_celltype_proportions(annotated)
  expect_s3_class(p2, "ggplot")

  p3 <- plot_marker_distributions(mif_cells, marker_map, cutoffs = mif_cutoffs)
  expect_s3_class(p3, "ggplot")
})
test_that("Phase 2 QC plots return ggplot objects", {
  skip_if_not_installed("ggplot2")
  data("mif_cells")
  data("mif_cutoffs")
  data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  fov <- assign_fov(annotated)
  qc <- qc_fov_summary(fov)

  expect_s3_class(plot_fov_distribution(fov), "ggplot")
  expect_s3_class(plot_qc_summary(qc), "ggplot")
  expect_s3_class(plot_metric_density(qc, metric = "DC_count"), "ggplot")
  p_dc <- plot_dc_count_distribution(qc)
  expect_s3_class(p_dc, "ggplot")
  expect_identical(p_dc$labels$x, "DC Count")
  expect_identical(p_dc$labels$y, "Density")
  expect_identical(p_dc$layers[[1]]$aes_params$fill, "steelblue")
  expect_equal(p_dc$layers[[1]]$aes_params$alpha, 0.5)
  expect_identical(p_dc$layers[[1]]$aes_params$colour, "black")
})
test_that("theme_mif returns a theme and works on a plot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(theme_mif(), "theme")
  p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3),
                       ggplot2::aes(x, y)) +
    ggplot2::geom_point() + theme_mif()
  expect_s3_class(p, "ggplot")
})

test_that("plot_cell_map / plot_fov_distribution support sample_id", {
  skip_if_not_installed("ggplot2")
  data("mif_cells"); data("mif_cutoffs"); data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  fov <- assign_fov(annotated)
  p1 <- plot_cell_map(annotated, sample_id = "C1001")
  expect_s3_class(p1, "ggplot")
  p2 <- plot_fov_distribution(fov, sample_id = "C1001")
  expect_s3_class(p2, "ggplot")
  expect_error(plot_cell_map(annotated, sample_id = "NOPE"))
})

test_that("plot_cell_maps_by_sample and plot_fov_distributions_by_sample write files", {
  skip_if_not_installed("ggplot2")
  data("mif_cells"); data("mif_cutoffs"); data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  fov <- assign_fov(annotated)
  d <- tempfile("mifmaps"); dir.create(d)
  p1 <- plot_cell_maps_by_sample(annotated, d)
  p2 <- plot_fov_distributions_by_sample(fov, d)
  expect_equal(length(p1), 3L)
  expect_equal(length(p2), 3L)
  expect_true(all(file.exists(c(p1, p2))))
  unlink(d, recursive = TRUE)
})

test_that("plot_group_comparison supports box and violin geoms", {
  skip_if_not_installed("ggplot2")
  df <- data.frame(
    Group = rep(c("R", "NR"), each = 10),
    value = c(rnorm(10, 10, 2), rnorm(10, 14, 2))
  )
  expect_s3_class(plot_group_comparison(df, "value", "Group", geom = "box"), "ggplot")
  expect_s3_class(plot_group_comparison(df, "value", "Group", geom = "violin"), "ggplot")
  expect_s3_class(plot_group_comparison(df, "value", "Group", geom = "both"), "ggplot")
  p <- plot_group_comparison(df, "value", "Group", geom = "both")
  expect_equal(length(p$layers), 3L)
  expect_identical(p$labels$y, "% of DC 20µm neighbors")
})

test_that("plot_celltype_radar draws an fmsb radar and excludes Other", {
  skip_if_not_installed("fmsb")
  data("mif_cells"); data("mif_cutoffs"); data("mif_rules"); data("mif_groups")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  annotated <- add_groups(annotated, mif_groups)

  f <- tempfile(fileext = ".png")
  df <- plot_celltype_radar(annotated, file = f, width = 8, height = 8, res = 100)
  expect_true(file.exists(f))
  expect_gt(file.info(f)$size, 0)
  # returned data frame: max/min rows + one row per group, columns = cell types
  expect_equal(rownames(df), c("max", "min", "NR", "R"))
  # "Other" removed by default
  expect_false("Other" %in% colnames(df))
  # percentages: each group row sums to 100
  s <- rowSums(df[c("NR", "R"), ])
  expect_true(all(abs(s - 100) < 1e-6))
  unlink(f)

  # exclude = NULL keeps Other
  f2 <- tempfile(fileext = ".png")
  df2 <- plot_celltype_radar(annotated, exclude = NULL, file = f2,
                             width = 8, height = 8, res = 100)
  expect_true("Other" %in% colnames(df2))
  unlink(f2)
})
test_that("plot_nearest_rank returns a ggplot", {
  skip_if_not_installed("ggplot2")
  data("mif_cells"); data("mif_cutoffs"); data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  fov <- assign_fov(annotated)
  rank <- nearest_neighbor_rank(find_neighbors(fov, radius = 60))
  p <- plot_nearest_rank(rank, top_n = 3)
  expect_s3_class(p, "ggplot")
  p2 <- plot_nearest_rank(rank, top_n = 3, by_sample = TRUE)
  expect_s3_class(p2, "ggplot")
  expect_error(plot_nearest_rank(rank, top_n = 0))
})

test_that("plot_neighbor_ring_radar accepts custom ring and group colours", {
  d <- expand.grid(
    Sample = paste0("S", 1:6),
    CellType = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
  d$mean_prop <- seq(0.1, 0.6, length.out = nrow(d))
  p <- plot_neighbor_ring_radar(
    d,
    group_colors = c(R = "#AA0000", NR = "#0000AA"),
    ring_colors = c(A = "#EEEEEE", B = "#DDDDDD", C = "#CCCCCC")
  )
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p$coordinates, "CoordPolar"))
  expect_true(inherits(p$coordinates, "CoordFixed"))
  expect_true(any(vapply(p$layers, function(x) {
    inherits(x$geom, "GeomPolygon")
  }, logical(1))))
  expect_error(
    plot_neighbor_ring_radar(d, group_colors = "red"),
    "one colour per group"
  )
})

test_that("plot_neighbor_comparison labels the selected method", {
  d <- expand.grid(
    Sample = paste0("S", 1:6), FOV = paste0("F", 1:3),
    CellType = c("cDC", "CD8Tn"), stringsAsFactors = FALSE
  )
  d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
  d$mean_prop <- seq_len(nrow(d)) / nrow(d)
  stat <- compare_neighbor_groups(
    d, "FOV", "cluster_permutation", group_order = c("R", "NR")
  )
  p <- plot_neighbor_comparison(d, stat, title = "test")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "sample-block permutation")
  expect_true(any(vapply(p$layers, function(x) {
    inherits(x$geom, "GeomText")
  }, logical(1))))
})
