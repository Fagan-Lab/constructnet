#' Base reconstruction
#'
#' @param series Numeric
#' @param ... List, optional parameters
#' @return igraph object
#' @export
fit <- function(series, ...) UseMethod("fit")

#' @export
fit.default <- function(series, ...) {
  stop("series not provided")
}

#' @export
fit.numeric <- function(series, ...) {
  options <- list(...)
  verbose <- length(options$verbose) && as.logical(options$verbose)

  if (verbose) {
    print("verbose")
  }

  igraph::random.graph.game(1, 1)
}
