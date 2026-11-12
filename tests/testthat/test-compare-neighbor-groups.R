make_grouped_fov_summary <- function() {
  d <- expand.grid(
    Sample = paste0("S", 1:6),
    FOV = paste0("F", 1:3),
    CellType = c("cDC", "CD8Tn"),
    stringsAsFactors = FALSE
  )
  d$Group <- ifelse(d$Sample %in% paste0("S", 1:3), "R", "NR")
  sample_effect <- ifelse(d$Group == "R", 0.2, 0)
  d$mean_prop <- 0.1 + sample_effect +
    as.numeric(factor(d$FOV)) * 0.005
  d
}

test_that("all three neighbour comparison methods are selectable", {
  d <- make_grouped_fov_summary()
  sample_result <- compare_neighbor_groups(
    d, "sample", "sample_wilcox", group_order = c("R", "NR")
  )
  fov_result <- compare_neighbor_groups(
    d, "FOV", "fov_wilcox", group_order = c("R", "NR")
  )
  cluster_result <- compare_neighbor_groups(
    d, "FOV", "cluster_permutation", group_order = c("R", "NR")
  )

  expect_true(all(sample_result$analysis_level == "sample"))
  expect_true(all(fov_result$method == "fov_wilcox"))
  expect_true(all(cluster_result$method == "cluster_permutation"))
  expect_true(all(cluster_result$permutation_unit == "biological sample"))
  expect_true(all(cluster_result$n_permutations == choose(6, 3)))
  expect_true(all(cluster_result$effect > 0))
  expect_true(all(c("p.value", "p.adjust", "method_label") %in%
                    names(cluster_result)))
})

test_that("FOV and sample methods report their actual unit counts", {
  d <- make_grouped_fov_summary()
  sample_result <- compare_neighbor_groups(d, "sample", "sample_wilcox")
  fov_result <- compare_neighbor_groups(d, "FOV", "fov_wilcox")

  expect_true(all(sample_result$n_sample1 == 3L))
  expect_true(all(sample_result$n_sample2 == 3L))
  expect_true(all(fov_result$n_fov1 == 9L))
  expect_true(all(fov_result$n_fov2 == 9L))
})

test_that("incompatible level and method combinations fail clearly", {
  d <- make_grouped_fov_summary()
  expect_error(
    compare_neighbor_groups(d, "sample", "fov_wilcox"),
    "requires"
  )
  expect_error(
    compare_neighbor_groups(d, "FOV", "sample_wilcox"),
    "requires"
  )
})
