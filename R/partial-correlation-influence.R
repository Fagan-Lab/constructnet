#' Reconstruction of graphs using the partial correlation influence.
#'
#'
#' @param TS Input matrix, X * L matrix consisting of L observations from N seneors.
#' @param index An index variable or set of index variables, which are assumed
#' to be confounders of all other variables.
#' @param threshold_type Which thresholding function to use on the matrix of weights
#' @param ... arguments
#'
#' @return A reconstructed graph with N nodes
#' @export
partial_correlation_influence <- function(TS, index=None, threshold_type='range', ...) {
  data <- t(TS)
  N = nrow(data)

  mask = matrix(1, N)
  if (!is.null(index)) {
    mask[index] <- FALSE
  }

  p_corr = matrix(NA, N, N)
  # FIX #
  p_corr = partial_corr(data[, mask], )
}