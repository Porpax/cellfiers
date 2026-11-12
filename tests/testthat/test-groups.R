test_that("add_groups merges group information", {
  data("mif_groups")
  d <- data.frame(sample_id = c("C1001", "C1002", "C9999"), v = 1:3)
  out <- add_groups(d, mif_groups)
  expect_equal(out$Group, c("R", "NR", NA))
})

test_that("compare_groups runs wilcox and t tests", {
  set.seed(1)
  df <- data.frame(
    Group = rep(c("R", "NR"), each = 10),
    value = c(rnorm(10, 10, 1), rnorm(10, 14, 1))
  )
  res <- compare_groups(df, value_col = "value", group_col = "Group")
  expect_true(all(c("group1", "group2", "statistic", "p.value") %in% names(res)))
  expect_true(res$p.value[1] < 0.05)

  res_t <- compare_groups(df, value_col = "value", group_col = "Group",
                          method = "t")
  expect_true(res_t$p.value[1] < 0.05)
})

test_that("compare_groups supports a by column", {
  set.seed(2)
  df <- data.frame(
    Group = rep(c("R", "NR"), each = 10),
    CellType = rep(c("DC", "CD8"), 10),
    value = c(rnorm(10, 10, 1), rnorm(10, 14, 1))
  )
  res <- compare_groups(df, value_col = "value", group_col = "Group",
                        by = "CellType")
  expect_equal(nrow(res), 2L)
  expect_setequal(res$CellType, c("DC", "CD8"))
})

test_that("compare_groups errors on invalid input", {
  expect_error(compare_groups(data.frame(a = 1), value_col = "a",
                              group_col = "b"), "group column")
})

test_that("plot_group_comparison returns a ggplot", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  df <- data.frame(
    Group = rep(c("R", "NR"), each = 30),
    CellType = rep(rep(c("DC", "CD8"), each = 15), 2),
    value = c(rnorm(15, 10, 2), rnorm(15, 12, 2),
              rnorm(15, 8, 2), rnorm(15, 14, 2))
  )
  p <- plot_group_comparison(df, value_col = "value", group_col = "Group",
                             by = "CellType")
  expect_s3_class(p, "ggplot")
  p2 <- plot_group_comparison(df, value_col = "value", group_col = "Group")
  expect_s3_class(p2, "ggplot")
})

test_that("analyze_neighbors returns the expected components", {
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
  res <- analyze_neighbors(fov, radius = 60,
                           center_types = c("cDC", "cDC(Mature/Active)"))
  expect_named(res, c("neighbors", "fov_summary", "nearest_rank"))
  expect_true(all(c("n_neighbors", "nearestNGB") %in% names(res$neighbors)))
  expect_true("mean_prop" %in% names(res$fov_summary))
  expect_true("Rank" %in% names(res$nearest_rank))
})