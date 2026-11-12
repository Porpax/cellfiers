test_that("read_qupath reads a file and adds sample_id and cell_id", {
  tmp <- tempfile(fileext = ".txt")
  df <- data.frame(
    Image = "img1",
    Name = "PathCellObject",
    Parent = "SLN_Tregion",
    `Centroid X µm` = c(1, 2),
    `Cell: Opal 780 mean` = c(500, 600),
    check.names = FALSE
  )
  write.table(df, tmp, sep = "\t", row.names = FALSE, quote = FALSE)

  out <- read_qupath(tmp, sample_id = "C9999")
  expect_s3_class(out, "tbl_df")
  expect_true(all(out$sample_id == "C9999"))
  expect_equal(out$cell_id, c("C9999_1", "C9999_2"))
  expect_true("Centroid X µm" %in% names(out))
  unlink(tmp)
})

test_that("read_qupath extracts sample id from file name", {
  tmp <- file.path(tempdir(), "Shanxi_SLN_C12345.txt")
  write.table(data.frame(a = 1), tmp, sep = "\t", row.names = FALSE)
  out <- read_qupath(tmp)
  expect_equal(unique(out$sample_id), "C12345")
  unlink(tmp)
})

test_that("read_qupath_batch combines samples", {
  d1 <- data.frame(v = 1, check.names = FALSE)
  d2 <- data.frame(v = 2, check.names = FALSE)
  f1 <- file.path(tempdir(), "A_C1001.txt")
  f2 <- file.path(tempdir(), "B_C1002.txt")
  write.table(d1, f1, sep = "\t", row.names = FALSE)
  write.table(d2, f2, sep = "\t", row.names = FALSE)

  combined <- read_qupath_batch(tempdir(), pattern = "C100[12][.]txt$")
  expect_equal(nrow(combined), 2L)
  expect_setequal(combined$sample_id, c("C1001", "C1002"))

  lst <- read_qupath_batch(tempdir(), pattern = "C100[12][.]txt$", combine = FALSE)
  expect_true(is.list(lst))
  expect_equal(names(lst), c("C1001", "C1002"))
  unlink(c(f1, f2))
})