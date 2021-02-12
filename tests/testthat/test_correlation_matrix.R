context("Testing Network Distance")

test_that("correlation_matrix", {
  mydata <- matrix(
    c(2, 4, 3, 1, 5, 7, 3, 4, 7),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )

  output <- correlation_matrix_fit(mydata, cutoffs = c(-1, 1))
})
