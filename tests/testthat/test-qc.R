make_annotated_fov <- function() {
  data("mif_cells")
  data("mif_cutoffs")
  data("mif_rules")
  marker_map <- c(
    CD8 = "Cell: Opal 780 mean", CD11c = "Cell: Opal 570 mean",
    CD69 = "Cell: Opal 690 mean", CD3 = "Cell: Opal 480 mean",
    PDL1 = "Cell: Opal 620 mean", CD45RA = "Cell: Opal 520 mean"
  )
  annotated <- run_annotation(mif_cells, mif_cutoffs, mif_rules, marker_map)
  assign_fov(annotated)
}

test_that("qc_fov_summary returns the expected columns", {
  fov <- make_annotated_fov()
  qc <- qc_fov_summary(fov)
  expect_true(all(c("Sample", "FOV", "TotalCells", "DC_count", "Tcell_count",
                    "DC_ratio", "Tcell_ratio", "MeanCells", "LowCellFlag",
                    "HighDCRatioFlag", "InitialQC") %in% names(qc)))
  expect_true(all(qc$LowCellFlag %in% c("OK", "LowCellNum")))
  expect_true(all(qc$HighDCRatioFlag %in% c("OK", "HighDCRatio")))
  # total cells across FOVs equals number of cells
  expect_equal(sum(qc$TotalCells), nrow(fov))
  # DC counts are non-negative
  expect_true(all(qc$DC_count >= 0))
})

test_that("qc_fov_summary applies an adjustable high-DC-ratio rule", {
  df <- data.frame(
    sample_id = rep("S1", 20),
    FOV = rep(c("F1", "F2"), each = 10),
    CellType = c(rep("cDC", 8), rep("CD8Tn", 2),
                 rep("cDC", 6), rep("CD8Tn", 4))
  )
  qc <- qc_fov_summary(df, max_dc_ratio = 0.70)
  expect_equal(qc$HighDCRatioFlag[qc$FOV == "F1"], "HighDCRatio")
  expect_equal(qc$HighDCRatioFlag[qc$FOV == "F2"], "OK")

  disabled <- qc_fov_summary(df, max_dc_ratio = NULL)
  expect_true(all(disabled$HighDCRatioFlag == "NotEvaluated"))
  expect_error(qc_fov_summary(df, max_dc_ratio = 1.1), "max_dc_ratio")
})

test_that("qc_fov_summary honours custom dc_types / tcell_types", {
  fov <- make_annotated_fov()
  qc <- qc_fov_summary(fov, dc_types = "cDC", tcell_types = "CD3")
  # no CD3 type exists -> Tcell_count should be 0
  expect_true(all(qc$Tcell_count == 0))
})

test_that("filter_qc_fovs keeps only passing FOVs", {
  # deterministic scenario: one FOV has far fewer cells than the sample mean
  df <- data.frame(
    sample_id = rep("S1", 51),
    FOV = c(rep("S1_TR1_FOV1", 20), rep("S1_TR1_FOV2", 30),
            rep("S1_TR1_FOV3", 1)),
    CellType = "CD3"
  )
  qc <- qc_fov_summary(df)
  expect_true("S1_TR1_FOV3" %in% qc$FOV[qc$LowCellFlag == "LowCellNum"])
  clean <- filter_qc_fovs(df, qc)
  expect_setequal(unique(clean$FOV), c("S1_TR1_FOV1", "S1_TR1_FOV2"))
  expect_lt(nrow(clean), nrow(df))
})

test_that("remove_extreme_fovs trims extreme metric values", {
  qc <- data.frame(
    FOV = paste0("F", 1:10),
    DC_count = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 100),
    TotalCells = 50
  )
  qc_f <- remove_extreme_fovs(qc, metric = "DC_count")
  expect_lt(nrow(qc_f), nrow(qc))
  q <- quantile(qc$DC_count, c(0.05, 0.95))
  expect_true(all(qc_f$DC_count >= q[1] & qc_f$DC_count <= q[2]))

  cells <- data.frame(FOV = rep(qc$FOV, each = 2), v = 1:20)
  both <- remove_extreme_fovs(qc, metric = "DC_count", data = cells)
  expect_true(is.list(both))
  expect_true(all(both$data$FOV %in% both$qc$FOV))
})

test_that("remove_outliers_iqr removes outliers globally and per group", {
  df <- data.frame(
    group = rep(c("A", "B"), each = 12),
    value = c(rnorm(12, 5, 1), c(rnorm(11, 50, 10), 500))
  )
  clean <- remove_outliers_iqr(df, value_col = "value")
  expect_lt(nrow(clean), nrow(df))

  clean_g <- remove_outliers_iqr(df, value_col = "value", group_col = "group")
  expect_lt(nrow(clean_g), nrow(df))
  # the extreme 500 belongs to group B and must be gone
  expect_false(500 %in% clean_g$value)
})

test_that("run_fov_qc returns consistent Phase 2 results", {
  fov_input <- make_annotated_fov()
  expect_warning(
    result <- run_fov_qc(fov_input, probs = c(0, 1)),
    "already exists"
  )
  expect_named(result, c(
    "data_fov", "qc", "qc_pass", "qc_final", "data_final",
    "metric_quantiles"
  ))
  expect_true(all(result$qc_pass$LowCellFlag == "OK"))
  expect_true(all(result$qc_pass$HighDCRatioFlag == "OK"))
  expect_true(all(result$qc$FinalQC %in% c("OK", "Exclude")))
  expect_setequal(unique(result$data_final$FOV), result$qc_final$FOV)
  expect_lte(nrow(result$data_final), nrow(result$data_fov))
})
