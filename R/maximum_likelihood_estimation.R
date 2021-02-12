
#' Reconstruction of graphs using maximum likelihood estimation

#' @param TS Matrix consisting of L observations from N sensors.
#'
#' @param rate rate term in maximum likelihood
#' @param stop_criterion if True, prevent overly-long runtimes
#' @param threshold_type Which thresholding function to use on the matrix.
#' @param ... Arguments
#'
#' @export
maximum_likelihood_estimation_fit <- function(TS, rate = 1.0, stop_criterion = T, threshold_type = "degree", ...) {
  # Infer inter-node coupling weights using maximum likelihood estimation
  # methods.
  # The results dictionary also stores the weight matrix as
  # `'weights_matrix'` and the thresholded version of the weight matrix
  # as `'thresholded_matrix'`.

  # Returns
  # -------
  #   G
  #      a reconstructed graph.
  N <- nrow(TS)
  L <- ncol(TS)
  rate <- rate / L
  W <- matrix(0, N, N)
  s1 <- TS[, -L]

  nloop <- 10000
  for (i0 in 1:N) {
    st1 <- matrix((TS[i0, -1]), nrow = 1)
    w <- matrix(0, 1, N)
    h <- matrix(0, 1, (L - 1))
    cost <- matrix(100, 1, nloop)

    for (iloop in 1:nloop) {
      dw <- t(s1 %*% t((st1 - tanh(h))))

      w <- w + (rate * dw)
      h <- t(t(s1) %*% t(w))
      cost[iloop] <- mean((st1 - tanh(h))**2)

      if (stop_criterion && iloop != 1 && (cost[iloop] >= cost[iloop - 1])) {
        break
      }
    }
    W[i0, ] <- w
  }
  W_thresh <- threshold(W, threshold_type, ...)

  G <- create_graph(W_thresh)

  structure(
    list(
      weights_matrix = W,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "MaximumLikelihoodEstimation"
  )
}
