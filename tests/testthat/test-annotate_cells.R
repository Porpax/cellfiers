test_that("annotate_cells assigns types from hand-made rules", {
  rules <- data.frame(
    CellType = c("DC", "Tcell"),
    CD11c = c("+", "-"),
    CD3 = c("-", "+"),
    stringsAsFactors = FALSE
  )
  cells <- data.frame(
    CD11c = c("+", "-", "+", NA),
    CD3 = c("-", "+", "+", "-"),
    stringsAsFactors = FALSE
  )
  out <- annotate_cells(cells, rules)
  expect_equal(out$CellType, c("DC", "Tcell", "Other", "Other"))
})

test_that("annotate_cells accepts 1/0 and TRUE/FALSE statuses", {
  rules <- data.frame(
    CellType = c("DC", "Tcell"),
    CD11c = c("+", "-"),
    CD3 = c("-", "+"),
    stringsAsFactors = FALSE
  )
  cells <- data.frame(
    CD11c = c(1, 0, TRUE, FALSE),
    CD3 = c(0, 1, FALSE, TRUE)
  )
  out <- annotate_cells(cells, rules)
  expect_equal(out$CellType, c("DC", "Tcell", "DC", "Tcell"))
})

test_that("annotate_cells respects rule order (first match wins)", {
  rules <- data.frame(
    CellType = c("Broad", "Specific"),
    A = c("+", "+"),
    B = c(NA, "+"),
    stringsAsFactors = FALSE
  )
  cells <- data.frame(A = c("+", "+"), B = c("+", "-"))
  out <- annotate_cells(cells, rules)
  # cell 1 matches both rules -> Broad (first); cell 2 matches only Broad
  expect_equal(out$CellType, c("Broad", "Broad"))
})

test_that("annotate_cells uses a custom celltype column name", {
  rules <- data.frame(CellType = "DC", CD11c = "+", CD3 = "-")
  cells <- data.frame(CD11c = "+", CD3 = "-")
  out <- annotate_cells(cells, rules, celltype_col = "Name")
  expect_equal(out$Name, "DC")
})
