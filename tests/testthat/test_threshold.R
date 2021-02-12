context("Testing Network Distance")

test_that("threshold", {
  m <- matrix(c(1, 4, 3, 1, 5, 7, 3, 4, 7), nrow = 3, ncol = 3, byrow = TRUE)
  cutoffs <- list(c(-1, 1), c(1, 3))
  testthat::expect_error(threshold(m, "custom"))
})
