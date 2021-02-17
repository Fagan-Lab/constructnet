#' Graph reconstruction algorithm.
#'
#' @param TS Matrix consisting of L observations from N variables.
#'
#' @param alpha Coefficient of penalization, higher values means more sparseness
#' @param max_iter Maximum number of iterations.
#' @param tol Stop the algorithm when the duality gap is below a certain threshold.
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param ... Arguments
#'
#' @export
graphical_lasso_fit <- function(TS, alpha = 0.01, max_iter = 100, tol = 0.0001, threshold_type = "degree", ...) {
  emp_cov <- stats::cov(t(TS))
  prec <- glasso::glasso(emp_cov, alpha, maxit = max_iter, thr = tol)$wi
  cov <- glasso::glasso(emp_cov, alpha, maxit = max_iter, thr = tol)$w

  W_thresh <- threshold(cov, threshold_type, ...)

  G <- create_graph(W_thresh)

  structure(
    list(
      weights_matrix = cov,
      percision_matrix = prec,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "GraphicalLasso"
  )
}
