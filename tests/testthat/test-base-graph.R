context("Testing Network Distance")

### Testing Base --------------------------------------------------------------
test_that("Test base works with series", {
  g <- fit(1)

  expect_true(igraph::is_igraph(g))
})

test_that("Test base errors with non-numeric data", {
  expect_error(fit("a"))
})
