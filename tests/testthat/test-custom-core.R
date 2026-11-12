test_that("region selection supports exact, fixed and regex modes", {
  x <- data.frame(Parent = c("Tumor_core", "Tumor_edge", "Stroma"))
  expect_equal(nrow(filter_region(x, regions = "Stroma")), 1)
  expect_equal(nrow(filter_region(x, regions = "Tumor", match = "fixed")), 2)
  expect_equal(nrow(filter_region(
    x, regions = "tumou?r", match = "regex", ignore_case = TRUE
  )), 2)
})

test_that("generic target neighbourhood uses inverse-degree weighting", {
  x <- data.frame(
    sample_id = "S1",
    CellType = c("Centre", "Centre", "TargetA", "TargetB", "Other"),
    FOV = "S1_FOV1",
    x = c(0, 10, 5, 0, 100),
    y = c(0, 0, 0, 5, 100)
  )
  z <- analyze_target_neighborhoods(
    x, center_types = "Centre", target_types = c("TargetA", "TargetB"),
    other_types = "Other", radius = 20, grid_size = 200
  )
  expect_equal(nrow(z$center_level), 2)
  expect_equal(z$settings$shared_weight, "inverse_degree")
  expect_equal(
    as.vector(tapply(z$edges$EdgeWeight, z$edges$TargetCellRow, sum)),
    rep(1, length(unique(z$edges$TargetCellRow)))
  )
  expect_true(all(c("Weighted_TargetA", "Weighted_TargetB") %in%
                    names(z$fov_level)))
})

test_that("radius above 30 um warns but remains available", {
  x <- data.frame(
    sample_id = "S1", CellType = c("Centre", "Target"), FOV = "F1",
    x = c(0, 5), y = c(0, 0)
  )
  expect_warning(
    analyze_target_neighborhoods(
      x, center_types = "Centre", target_types = "Target", radius = 31
    ),
    "greater than 30"
  )
})
