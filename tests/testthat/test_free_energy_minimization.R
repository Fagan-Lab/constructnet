context("Testing Network Distance")

test_that("free_energy_minimization", {
  TS <- matrix(
    c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9),
    nrow = 3,
    ncol = 7,
    byrow = TRUE
  )

  f <- free_energy_minimization_fit(TS, avg_k = 1)
})
