#' Creates a minimum spanning tree connecting the sensors.
#'
#' @param inputArray inputArray
#' @param distance distance
#' @param ... ...
#'
#' @export
construct_correlation_spanning <- function(inputArray, distance, ...) {
  UseMethod("construct_correlation_spanning")
}

#' @export
construct_correlation_spanning.igraph <- function(inputArray, distance, ...) {
}
