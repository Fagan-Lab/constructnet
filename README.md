
# constructnet

<!-- badges: start -->
[![R build status](https://github.com/travisbyrum/constructnet/workflows/R-CMD-check/badge.svg)](https://github.com/travisbyrum/constructnet/actions)
<!-- badges: end -->

constructnet: Network reconstruction from time series.

## Installation

You can install the released version of constructnet from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("constructnet")
```

## Development

### Document

To document the package use:

``` r
devtools::document()
```

This will add any new function `.Rd` files and update existing functions.

### Check

To check that the package installs correctly and all tests are passing, use:

``` r
devtools::check(document = FALSE)
```

Alternatively, you can use the check button on the `Build` tab in Rstudio.

### Adding a dependency

To add a package dependency use:

``` r
usethis::use_package("package")
```

This will correctly update the `DESCRIPTION` with the new dependency.

### Adding a test

To add a test for a given function, add a corresponding test file in the `tests/testthat` directory.  To run all the tests in the package, use:

``` r
devtools::test()
```

This is also available as the `CMD + SHIFT + t` command in Rstudio.

### Install local package

To install package locally for development use:

``` r
devtools::load_all(".")
```

This is also available as the `CMD + SHIFT + l` command in Rstudio.
