
#' Graph reconstruction algorithm from time series data


#' @param TS input matrix, N * L matrix consisting of L observations
#'           from N sensors.
#'
#' @param tau Number of time steps for a single time-lag.
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param cutoffs cutoffs for threshold function
#' @param ... arguments
#'
#' @export
convergent_cross_mapping_fit <-function(TS, tau = 1, threshold_type='range', cutoffs=list(c(0.95, Inf)), ...) {
    #  Convergent cross-mapping infers dynamical causal relation between
    #  vairiables from time series data

    # Returns
    # -------
    #
    #   G
    # A reconstructed graph with :`N` nodes.

  data <- t(TS)
  L = nrow(data)
  N = ncol(data)

  if(L < 3 + (N - 1) * (1 +  tau)) {
    stop('Need more data. L must be not less than 3+(N-1)*(1+tau).')
  }


  # Create shadow data cloud for each variable
  shadows <- matrix(list(matrix(list(), nrow = 1, ncol = L)), nrow = 1, ncol = N)
  for (j in 1 : N) {
    shadows[[1, j]] = shadow_data_cloud(data[, j], N, tau)
  }

  # Obtain nearest neighbors of points in the shadow data clould and
  # their weights for time series estimates
  neighbors <- matrix(list(matrix(list(), nrow = 1, ncol = L)), nrow = 1, ncol = N)
  for(i in 1 : N) {
    neighbors[[1, i]] = nearest_neighbors(shadows[[i]], L)[[1]]
  }
  distances <- matrix(list(matrix(list(), nrow = 1, ncol = L)), nrow = 1, ncol = N)
  for(i in 1 : N) {
    distances[[1, i]] = nearest_neighbors(shadows[[i]], L)[[2]]
  }
  nei_weights <- matrix(list(matrix(list(), nrow = 1, ncol = L)), nrow = 1, ncol = N)
  for (i in 1 : N) {
    nei_weights[[1, i]] <- neighbor_weights(distances[[i]])
  }

  # For every variable X and every other variable Y,
  # construct the estimates of Y from X's shadow data cloud and
  # compute the Pearson correlation between Y and its estimates
  # along with the p-value
  correlation = matrix(1, N, N)
  pvalue = matrix(0, N, N)
  p <- gtools::permutations(N,2, 1:N)
  for(x in 1:nrow(p)){
    i <- p[x, ][1]
    j <- p[x, ][2]
    estimates <- time_series_estimates(data[, j], neighbors[[i]], nei_weights[[i]])
    M <- length(estimates)
    correlation[i, j] = cor.test(estimates, data[(L-M+1):L, j], method = "pearson")[[4]][[1]]
    pvalue[i, j] = cor.test(estimates, data[(L-M+1):L, j], method = "pearson")[[3]]
  }

  # Build the reconstructed graph by finding significantly correlated
  # variables,
  # where weights of reconstructed edges are selected such that we can
  # sort zero p-values in decreasing order and tell edges with zero p-value
  weights = pvalue
  weights[weights > 0] <- -log10(weights[weights > 0])
  weights[weights == 0] <- Inf
  A = threshold(weights, threshold_type, cutoffs = cutoffs, ...)
  G = create_graph(A)

  structure(
    list(
      correlation_matrix = correlation,
      pvalues_matrix = pvalue,
      weights_matrix = weights,
      thresholded_matrix = A,
      graph = G
    ),
    class = "ConvergentCrossMapping"
  )
}

shadow_data_cloud <- function(data, N, tau) {
    # Return the lagged-vector data cloud of a given variable's time series.
    #
    # Parameters
    # ----------
    # data : Length:`L` 1D array of a single variable's times series.
    #
    # N (int): Number of variables.
    #
    # tau (int): Number of time steps for a single time-lag.
    #
    # Returns
    # -------
    # shadow : :`M \times aN` array of the lagged-vector data cloud,
  L <- length(data)
  M = L - (N - 1) * tau     # Number of points in the shadow data cloud
  shadow <- matrix(0, M, N)

  for (j in N : 1) {            # Fill in column values from the right
    delta <- (N - j)*tau        # Amount of time-lag for this column
    shadow[, j] = data[(delta +  1) : (delta + M )]   # + 1 because of the difference between Python and R
  }

  results <- shadow
}


nearest_neighbors <- function(shadow, L) {
  # Return time indices of the N+1 nearest neighbors for every point in the
  # shadow data cloud and their corresponding Euclidean distances.
  #
  # Parameters
  # ----------
  #   shadow : Array of the shadow data cloud.
  #
  # L (int): Number of observations in the time series.
  #
  # Returns
  # -------
  #   nei : `M \times (N+1)` array of time indices of nearest
  # neighbors where :`M` is the number of points in the
  # shadow data cloud.
  #
  # dist : `M \times (N+1)` array of corresponding Euclidean
  # distance between the data point to the neighbors.

  M = nrow(shadow)
  N = ncol(shadow)
  K <- N + 1

  if(K < M / 2) {
    method = 'kd_tree'
  } else {
    method = 'brute'
  }

  nbrs = FNN::get.knn(shadow, k = K, algorithm = method)
  nei = nbrs[[1]]
  dist = nbrs[[2]]

  # # Remove the first column for both arrays
  # nei <- nei[, -c(1)]
  # dist <- dist[, -c(1)]

  nei = nei + L - M +  1
  mylist <- list("nei" = nei, "dist" = dist)
  return(mylist)
}


neighbor_weights <- function(dist){
    # Return the weights of neighbors in time seires estimates.
    #
    # Parameters
    # ----------
    # dist :`M \times (N+1)` array of Euclidean distances between a
    #                    point to its nearest neighbors in the shadow data cloud
    #                    (sorted by increasing order of distances), where `M` is
    #                    the number of points in the shadow data cloud.
    #
    # Returns
    # -------
    # wei:`M \times (N+1)` array of exponentially decaying weights
    #                   of the nearest neighbors.

  r <- nrow(dist)
  c <- ncol(dist)

  # For those where the nearest distance is positive, assign weights that are
  # exponentially decaying with the ratio to the nearest distance
  dist_pnd <- matrix(dist[dist[, 1] > 0], ncol = c )
  dist[(dist[, 1] > 0)] <- exp(dist[(dist[, 1] > 0)] / -dist_pnd[, 1])
  dist[(dist[, 1] < 0)] <- 0
  for(i in 1: r){
    if(dist[i][1] == 0){
      dist[i][1] <- 1
    }
  }
  dist <- dist / rowSums(dist)
  results <- dist
}


time_series_estimates <- function(data_y, nei_x, wei_x){
    # Return estimates of variable :math:`Y` from variable :math:`X`'s shadow data cloud.
    #
    # Parameters
    # ----------
    # data_y : 1D array of variable :math:`Y`'s time series data.
    #
    # nei_x : :math:`M \times (N+1)` array of time indices of nearest
    #                     neighbors in :math:`X`'s shadow data cloud, where :math:`M` is the
    #                     number of points in the shadow data cloud.
    #
    # wei_x : Array of corresponding weights of the nearest
    #                     neighbors in :math:`X`'s shadow data cloud.
    #
    # Returns
    # -------
    # ests : Length-:math:`M` 1D array of estimates of :math:`Y`'s time series.

  r = nrow(nei_x)
  c = ncol(nei_x)

  nei_x = nei_x - 1
  ests <- rowSums((matrix(data_y[nei_x], nrow = r, ncol = c)) * wei_x)

  ests
}
