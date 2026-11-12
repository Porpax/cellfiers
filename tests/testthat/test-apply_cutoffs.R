test_that("apply_cutoffs adds +/- marker columns per sample", {
  data("mif_cells")
  data("mif_cutoffs")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  out <- apply_cutoffs(mif_cells, mif_cutoffs, marker_map)
  expect_true(all(c("CD8", "CD11c", "CD69", "CD3", "PDL1", "CD45RA") %in% names(out)))
  expect_true(all(out$CD3 %in% c("+", "-")))
  # a cell with intensity >= cutoff must be "+"
  s <- "C1001"
  sub <- out[out$sample_id == s, ]
  cut <- mif_cutoffs[mif_cutoffs$SampleID == s, ]
  expect_true(all(sub$CD3[sub$`Cell: Opal 480 mean` >= cut$CD3] == "+"))
  expect_true(all(sub$CD3[sub$`Cell: Opal 480 mean` < cut$CD3] == "-"))
})

test_that("apply_cutoffs warns on missing sample cutoff", {
  data("mif_cells")
  data("mif_cutoffs")
  marker_map <- c(CD3 = "Cell: Opal 480 mean")
  cut <- mif_cutoffs[mif_cutoffs$SampleID == "C1001", ]
  expect_warning(apply_cutoffs(mif_cells, cut, marker_map), "No cutoff row")
})

test_that("apply_cutoffs handles custom pos/neg labels", {
  data("mif_cells")
  data("mif_cutoffs")
  marker_map <- c(CD3 = "Cell: Opal 480 mean")
  out <- apply_cutoffs(mif_cells, mif_cutoffs, marker_map, pos = "pos", neg = "neg")
  expect_true(all(out$CD3 %in% c("pos", "neg")))
})