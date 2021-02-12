
#' Graph reconstruction algorithm.


#' @param TS N * L matrix consisting of L observations from N sensors.
#'
#' @param remove_largest If False, all the eigenvectors associated to the
#' significant eigenvalues will be used to reconstruct the
#' de-noised empirical correlation matrix. If True, the
#' eigenvector associated to the largest eigenvalue is going to be excluded from
#' the recontruction step.
#' @param metric_distance If False, a signed graph is obtained. If True, the correlation
#' is transformed by defining a metric distance between each pair of nodes
#' @param tol avoid minor difference between python and r in case of calculating eigen values
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param ... Arguments
#'
#' @export
marchenko_pastur_fit <- function(TS, remove_largest = F, metric_distance = F, tol = 1e-15, threshold_type = "range", ...) {
  # Create a correlation-based graph using Marchenko-Pastur law to remove noise.

  # Returns
  # -------
  #
  #   G
  #     A reconstructed graph.
  N <- nrow(TS)
  L <- ncol(TS)
  if (N > L) {
    stop('"L must be greater or equal than N."')
  }

  Q <- L / N
  C <- cor(t(TS))
  ev <- eigen(t(C))
  w <- ev$values
  v <- ev$vectors
  w_min <- 1 + 1 / Q - 2 * sqrt(1 / Q)
  w_max <- 1 + 1 / Q + 2 * sqrt(1 / Q)

  selected <- (w < w_min + tol) | (w > w_max - tol)
  if (sum(selected) == 0) {
    G <- igraph::make_empty_graph(n = N)
    result <- G
    return(result)
  }

  if (remove_largest) {
    selected[1] <- F
  }

  w_signal <- as.matrix(w[selected])
  v_signal <- as.matrix(v[, selected])
  C_signal <- v_signal %*% (diag(w_signal)) %*% t(v_signal)

  if (metric_distance) {
    C_signal <- sqrt(2 * (1 - C_signal))
  }

  t <- threshold(C_signal, threshold_type, ...)

  G <- create_graph(t)

  structure(
    list(
      weights_matrix = C_signal,
      thresholded_matrix = t,
      graph = G
    ),
    class = "GraphicalLasso"
  )
}
