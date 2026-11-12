test_that("filter_region keeps and drops matching rows", {
  data("mif_cells")
  tzone <- filter_region(mif_cells, regions = "SLN_Tregion")
  expect_true(all(tzone$Parent == "SLN_Tregion"))
  expect_lt(nrow(tzone), nrow(mif_cells))

  outside <- filter_region(mif_cells, regions = "SLN_Tregion", keep = FALSE)
  expect_true(all(outside$Parent != "SLN_Tregion"))
  expect_equal(nrow(tzone) + nrow(outside), nrow(mif_cells))
})

test_that("filter_region errors when parent column missing", {
  expect_error(filter_region(data.frame(a = 1), regions = "x"))
})