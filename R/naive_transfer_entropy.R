<<<<<<< HEAD
#' Implementation for graph reconstruction algorithm.
#'
#'
#' @param TS The input vector with 'N' sensors
#' @param delay_max The number of timesteps in the past to aggregate and average
#' @param n_bins The number of bins to turn values in the time series to categorical
#' data, which is a pre-processing step to compute entropy.
#' @param threshold_type Specify which thresholding function to use on the matrix
#' of weights.
#' @param ... arguments
#' @return A reconstructed graph with 'N' nodes.
#'
#' @export
naive_transfer_entropy <- function(TS, delay_max, n_bins, threshold_type, ...) {
  # The resulting network is asymmetric, and each element represents the amount of
  # information contained about the future states by knowing the past states.

  N = ncol(TS)
  L = nrow(TS)
  data <- t(TS)

  if(delay_max >= L) {
    stop("Max steps of delay exceeds time series length.")
  }

  data <- categorized_data(data, n_bins)
  # Compute the transfer entropy of every tuple of nodes
  TE = matrix(0, N, N)

  p <- gtools::permutations(N, 2, 1:N)
  for(x in 1:nrow(p)) {
    # Check several delay values and average them together
    # This average is naive, but appears to be sufficient in
    # some circumstances
    i <- p[x, ][1]
    j <- p[x, ][2]
    for(delay in 1: delay_max + 1) {
      te_list = transfer_entropy(data[, i], data[, j], delay)
    }
    TE[i, j] = mean(te_list)
  }

  TE_thresh = threshold(TE, threshold_type, ...)
  G = create_graph(TE_thresh)

  structure(
    list(
      weights.matrix = TE
      thresholded_matrix = TE_tresh
=======
#' Graph reconstruction algorithm
#' Calculates the transfer entropy from i --> j.
#' The resulting network is asymmetric
#'
#' @param TS  array consisting of L observations from N sensors.
#' @param delay_max the number of timesteps in the past to aggregate and average in order to get TE_{ij}
#' @param n_bins the number of bins to turn values in the time series to categorical
#'  data, which is a pre-processing step to compute entropy.
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param ... arguments
#'
#' @return a reconstructed graph with N nodes.
#' @export
#'
naive_transfer_entropy_fit <- function(TS, delay_max=1, n_bins=2, threshold_type='range', ...) {
  
  N = nrow(TS)
  L = ncol(TS)
  data = t(TS)
  if(delay_max >= L) {
    stop("Max steps of delay exceeds time series length.")
  }
  
  data = categorized_data(data, n_bins)
  TE = matrix(0, N, N)
  p <- gtools::permutations(N, 2, 1:N)
  te_list = vector()
  for (x in 1:nrow(p)) {
    i <- p[x, ][1]
    j <- p[x, ][2]
    for (delay in 1:(delay_max )) {
      te_list = append(te_list, transfer_entropy(data[,i], data[,j], delay))
    }
    TE[i, j] = mean(te_list)
    te_list = vector()
  }
  
  TE_thresh = threshold(TE, threshold_type, ...)
  G = create_graph(TE_thresh)
  
  structure(
    list(
      weights_matrix = TE,
      thresholded_matrix = TE_thresh,
>>>>>>> d12faab50d7ea1a6c0a6844cd486fca7746138c7
      graph = G
    ),
    class = "NaiveTransferEntropy"
  )
<<<<<<< HEAD
}

transfer_entropy <- function(X, Y, delay) {
=======
  
}



transfer_entropy <- function(X, Y, delay){
>>>>>>> d12faab50d7ea1a6c0a6844cd486fca7746138c7
  # This is a TE implementation: asymmetric statistic measuring the reduction
  # in uncertainty for the dynamics of Y given the history of X. Or the
  # amount of information from X to Y. The calculation is done via conditional
  # mutual information.
<<<<<<< HEAD
  #
  # Parameters
  #----------
  # X (np.ndarray): time series of categorical values from node :math:`i`
  # Y (np.ndarray): time series of categorical values from node :math:`j`
  # delay (int): steps with which node :math:`i` past state is accounted
  # Returns
  # -------
  # te (float): the transfer entropy from nodes i to j
  #
  X_past <- X[1:(length(X)-delay)]
  X_past <- as.matrix(X_past)

  Y_past <- Y[1:(length(Y)-delay)]
  Y_past <- as.matrix(Y_past)

  joint_past <- cbind(X_past, Y_past)

  Y_future <- Y[delay:length(Y)]
  Y_future <- as.matrix(Y_future)

  te = conditional_entropy(Y_future, Y_past)
  te = te - conditional_entropy(Y_future, joint_past)

  te
}
=======
  
  # Parameters
  # ----------
  # X : time series of categorical values from node i
  # Y : time series of categorical values from node j
  # delay (int): steps with which node i past state is accounted
  #
  # Returns
  # -------
  #   te (float): the transfer entropy from nodes i to j
  
  xr = length(X)
  yr = length(Y)
  X_p = t(t(as.matrix(X)[-((xr - delay + 1):xr),]))
  Y_p = t(t(as.matrix(Y)[-((yr - delay + 1):yr),]))
  joint = cbind(Y_p, X_p)
  
  Y_f = t(t(as.matrix(Y)[-(1:delay),]))
  te = conditional_entropy(Y_f, Y_p)
  te = te - conditional_entropy(Y_f, joint)
  
  te
}
>>>>>>> d12faab50d7ea1a6c0a6844cd486fca7746138c7
