context("Testing Network Distance")

test_that("marchenko_pastur", {
  TS <- matrix(
    c(6, 3, 1, 5, 3, 0, 5, 1, 1, 5, 6, 2, 5, 1, 2, 2),
    nrow = 4,
    ncol = 4
  )

  m <- marchenko_pastur_fit(
    TS,
    remove_largest = TRUE,
    metric_distance = TRUE,
    tol = 1e-15,
    threshold_type = "range",
    cutoffs = list(c(-1, 1))
  )
})
