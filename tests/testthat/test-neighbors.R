test_that("find_neighbors counts neighbours within radius", {
  d <- data.frame(
    sample_id = c("S1", "S1", "S1", "S1", "S1"),
    FOV = c("F1", "F1", "F1", "F2", "F2"),
    CellType = c("A", "B", "A", "A", "A"),
    x = c(0, 0, 100, 0, 0),
    y = c(0, 5, 100, 0, 10),
    cell_id = paste0("c", 1:5)
  )
  nb <- find_neighbors(d, radius = 20)

  # F1 cell1 (0,0) type A: exactly one neighbour (B at distance 5)
  expect_equal(nb$n_neighbors[1], 1L)
  expect_equal(nb$B[1], 1L)
  expect_equal(nb$A[1], 0L)
  expect_equal(nb$nearestNGB[1], "B")

  # F1 cell3 (100,100): isolated -> no neighbours
  expect_equal(nb$n_neighbors[3], 0L)
  expect_equal(nb$nearestNGB[3], "None")

  # F2: two A cells 10 um apart
  expect_equal(nb$n_neighbors[4], 1L)
  expect_equal(nb$A[4], 1L)
  expect_equal(nb$nearestNGB[4], "A")

  # row order preserved
  expect_identical(nb$cell_id, d$cell_id)
})

test_that("find_neighbors can include the centre cell", {
  d <- data.frame(
    sample_id = "S1", FOV = "F1",
    CellType = c("A", "B"),
    x = c(0, 10), y = c(0, 0),
    cell_id = c("c1", "c2")
  )
  nb <- find_neighbors(d, radius = 20, include_self = TRUE)
  expect_equal(nb$n_neighbors, c(2L, 2L))
  expect_equal(nb$A[1], 1L)  # self
  expect_equal(nb$B[1], 1L)
})

test_that("distinct cells at identical coordinates remain neighbours", {
  d <- data.frame(
    sample_id = "S1", FOV = "F1", CellType = c("A", "B"),
    x = c(0, 0), y = c(0, 0), cell_id = c("c1", "c2")
  )
  nb <- find_neighbors(d, radius = 20)
  expect_equal(nb$n_neighbors, c(1L, 1L))
  expect_equal(nb$nearestNGB, c("B", "A"))
})

test_that("neighbor_proportions computes per-centre proportions", {
  d <- data.frame(
    sample_id = c("S1", "S1", "S1"),
    FOV = c("F1", "F1", "F1"),
    CellType = c("DC", "DC", "T"),
    x = c(0, 0, 5),
    y = c(0, 5, 0),
    cell_id = paste0("c", 1:3)
  )
  nb <- find_neighbors(d, radius = 20)
  prop <- neighbor_proportions(nb, center_types = "DC")
  expect_equal(nrow(prop), 2L)
  expect_equal(prop$total_neighbors, c(2, 2))
  expect_equal(prop$Pro_DC, c(0.5, 0.5))
  expect_equal(prop$Pro_T, c(0.5, 0.5))
})

test_that("neighbor_fov_summary aggregates per FOV (long and wide)", {
  d <- data.frame(
    sample_id = c("S1", "S1", "S1", "S2", "S2"),
    FOV = c("F1", "F1", "F1", "F2", "F2"),
    CellType = c("DC", "DC", "T", "DC", "T"),
    x = c(0, 0, 5, 0, 10),
    y = c(0, 5, 0, 0, 0),
    cell_id = paste0("c", 1:5)
  )
  nb <- find_neighbors(d, radius = 20)
  summ <- neighbor_fov_summary(nb, center_types = "DC")
  expect_true(all(c("Sample", "FOV", "CellType", "mean_count", "mean_prop",
                    "n_center") %in% names(summ)))
  expect_equal(unique(summ$FOV), c("F1", "F2"))
  expect_equal(summ$n_center[summ$FOV == "F1" & summ$CellType == "DC"], 2L)
  wide <- neighbor_fov_summary(nb, center_types = "DC", wide = TRUE)
  expect_equal(nrow(wide), 2L)
  expect_true("mean_prop_DC" %in% names(wide))

  sample_summ <- neighbor_sample_summary(summ)
  expect_equal(nrow(sample_summ), 4L)
  expect_true(all(c("n_fov", "n_center") %in% names(sample_summ)))
})

test_that("nearest_neighbor_rank ranks by count per FOV", {
  d <- data.frame(
    sample_id = c("S1", "S1", "S1", "S1"),
    FOV = c("F1", "F1", "F1", "F1"),
    CellType = c("A", "B", "A", "B"),
    x = c(0, 3, 6, 9),
    y = c(0, 0, 0, 0),
    cell_id = paste0("c", 1:4)
  )
  nb <- find_neighbors(d, radius = 20)
  rank <- nearest_neighbor_rank(nb)
  expect_equal(nrow(rank), 2L)
  expect_equal(sum(rank$Count), 4L)
  expect_equal(sum(rank$Proportion), 1)
  expect_true(all(rank$Rank == 1))

  centers <- nearest_neighbor_rank(nb, center_types = "A")
  expect_equal(sum(centers$Count), 2L)
})
