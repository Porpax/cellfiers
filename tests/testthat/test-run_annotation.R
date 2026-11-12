test_that("run_annotation performs the full Phase-1 pipeline", {
  data("mif_cells")
  data("mif_cutoffs")
  data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  out <- run_annotation(
    mif_cells, mif_cutoffs, mif_rules, marker_map,
    regions = "SLN_Tregion"
  )
  expect_true("CellType" %in% names(out))
  expect_true(all(out$Parent == "SLN_Tregion"))
  valid <- c(mif_rules$CellType, "Other")
  expect_true(all(out$CellType %in% valid))
  # CD4Teff(Resting) must have CD3 = "+" and CD8 = "-"
  cd4r <- out[out$CellType == "CD4Teff(Resting)", ]
  expect_true(all(cd4r$CD3 == "+"))
  expect_true(all(cd4r$CD8 == "-"))
})

test_that("summarize_celltypes totals match the data", {
  data("mif_cells")
  data("mif_cutoffs")
  data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  out <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  summ <- summarize_celltypes(out)
  expect_equal(nrow(summ), 3L)
  expect_equal(sum(summ$Total), nrow(out))
  expect_true(all(c("Sample", "Total") %in% names(summ)))
  expect_true(all(mif_rules$CellType %in% names(summ)))
})