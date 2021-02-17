#' Graph reconstruction algorithm
#'
#' @param TS Array consisting of L observations from N variables.
#'
#' @param lag Time lag to consider.
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param ... Arguments
#'
#' @export
granger_causality_fit <- function(TS, lag = 1, threshold_type = "range", ...) {
  n <- nrow(TS)
  W <- matrix(0, n, n)

  for (i in 1:n) {
    xi <- split_data(TS[i, ], lag)$inputs
    yi <- split_data(TS[i, ], lag)$targets

    for (j in 1:n) {
      xj <- split_data(TS[j, ], lag)$inputs
      xij <- cbind(xi, xj)
      reg1 <- stats::lm(as.vector(yi) ~ xi)
      reg2 <- stats::lm(as.vector(yi) ~ xij)
      err1 <- yi - stats::predict(reg1, as.data.frame(xi))
      err2 <- suppressWarnings(yi - stats::predict(reg2, as.data.frame(xij)))

      std_i <- pracma::std(err1, 1)
      std_ij <- pracma::std(err2, 1)

      if (std_i == 0) {
        W[j, i] <- -99999999
      } else if (std_ij == 0) {
        W[j, i] <- 99999999
      } else {
        W[j, i] <- log(std_i) - log(std_ij)
      }
    }
  }

  # threshold the network
  W_thresh <- threshold(W, threshold_type, ...)

  # construct the network
  G <- create_graph(W_thresh)

  structure(
    list(
      weights_matrix = W,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "GrangerCausality"
  )
}


split_data <- function(TS, lag) {
  T <- length(TS)
  inputs <- matrix(0, (T - lag - 1), lag)
  targets <- matrix(0, (T - lag - 1))

  for (t in 1:((T - lag - 1))) {
    inputs[t, ] <- TS[lag + t - 1]
    targets[t] <- TS[t + lag]
  }

  list("inputs" = inputs, "targets" = t(targets))
}
