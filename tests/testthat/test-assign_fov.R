test_that("assign_fov adds FOV labels and preserves row order", {
  data("mif_cells")
  data("mif_cutoffs")
  data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  out <- assign_fov(annotated)
  expect_true("FOV" %in% names(out))
  expect_true(all(grepl("^C[0-9]+_TR[0-9]+_FOV[0-9]+$", out$FOV)))
  # row order preserved
  expect_identical(out$cell_id, annotated$cell_id)
  # two well-separated clusters -> two tissue regions per sample
  trs <- gsub(".*_TR([0-9]+)_FOV[0-9]+$", "\\1", out$FOV[out$sample_id == "C1001"])
  expect_equal(length(unique(trs)), 2L)
})

test_that("assign_fov works with dot-named centroid columns", {
  d <- data.frame(
    sample_id = "S1",
    `Centroid X µm` = c(50, 150, 950, 1050),
    `Centroid Y µm` = c(50, 150, 950, 1050),
    check.names = FALSE
  )
  out <- assign_fov(d, grid_size = 200, min_distance = 500)
  expect_true(all(grepl("_TR[0-9]+_FOV[0-9]+$", out$FOV)))
})

test_that("assign_fov handles a single grid (one region)", {
  d <- data.frame(
    sample_id = "S1",
    `Centroid X µm` = c(10, 20, 30),
    `Centroid Y µm` = c(10, 20, 30),
    check.names = FALSE
  )
  out <- assign_fov(d, grid_size = 200, min_distance = 500)
  expect_true(all(grepl("_TR1_FOV1$", out$FOV)))
})

test_that("assign_fov rejects unusable inputs", {
  empty <- data.frame(sample_id = character(), x = numeric(), y = numeric())
  expect_error(assign_fov(empty, x = "x", y = "y"), "at least one")

  bad <- data.frame(sample_id = "S1", x = NA_real_, y = 1)
  expect_error(assign_fov(bad, x = "x", y = "y"), "finite numeric")
})
