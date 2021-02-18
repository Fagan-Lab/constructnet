#' Reconstruction of graphs by minimizing a free energy of your data
#'
#' @param TS Matrix consisting of L observations from N different variables.
#'
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param ... Arguments
#'
#' @export
free_energy_minimization_fit <- function(TS, threshold_type = "degree", ...) {
  N <- nrow(TS)
  L <- ncol(TS)
  m <- rowMeans(TS[, -L]) # model average
  ds <- t(TS[, -L])
  ds <- apply(ds, 1, function(x) x - m)
  ds <- t(ds) # discrepancy
  t1 <- L - 1 # time limit

  # covariance of the discrepeancy
  c <- stats::cov.wt(ds, method = "ML")$cov
  c_inv <- solve(c)
  dst <- t(ds) # discrepancy at time t
  W <- matrix(1, N, N)
  nloop <- 10000 # failsafe

  for (i0 in 1:N) {
    TS1 <- TS[i0, -1]
    h <- TS1
    cost <- matrix(100, 1, nloop)

    for (iloop in seq_len(nloop)) {
      h_av <- mean(h) # average local field
      hs_av <- t(dst %*% (h - h_av) / t1) # deltaE_i delta\sigma_k
      w <- (hs_av %*% c_inv) # expectation under model
      h <- as.vector(t(t(TS[, -L]) %*% w[, ])) # estimate of local field
      TS_model <- tanh(h) # under kinetic Ising model
      # discrepancy cost
      cost[iloop] <- mean((TS1 - TS_model)**2)

      if (iloop != 1 && (cost[iloop] >= cost[iloop - 1])) {
        break
      }

      h <- h * (TS1 / TS_model)
    }

    W[i0, ] <- w[, ]
  }

  # threshold the network
  W_thresh <- threshold(W, threshold_type, ...)

  G <- create_graph(W_thresh)
  # construct the network
  structure(
    list(
      weights_matrix = W,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "FreeEnergyMinimization"
  )
}
