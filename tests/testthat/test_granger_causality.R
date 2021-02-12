context("Testing Network Distance")

test_that("granger_causality", {
  TS <- matrix(
    c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9),
    nrow = 3,
    ncol = 7,
    byrow = TRUE
  )

  s <- granger_causality_fit(TS, cutoffs = c(-1, 1))
})
