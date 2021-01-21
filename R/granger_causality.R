# granger_causality.R
# --------------
# status: Finished draft and some simple tests
# Graph reconstruction algorithm
# author: Charles Murphy
# converted by: Zhaoyi Zhuang

# source(here::here('constructnet', 'R', 'threshold.R'))
# source(here::here('constructnet', 'R', 'graph.R'))

granger_causality_fit <- function(TS, lag = 1, threshold_type = 'range', ...){
    # Reconstruct a network based on the Granger causality.
    # It reconstructs the network by calculating the Granger
    # causality for each pair of nodes.
    #
    # Parameters
    # ----------
    #
    #   TS matrix
    # Array consisting of :`L` observations from :`N`
    # sensors.
    #
    # lag (int)
    # Time lag to consider.
    #
    # threshold_type (str)
    # Which thresholding function to use on the matrix of
    # weights.
    #
    # Returns
    # --------
    #
    #   G
    # A reconstructed graph with :`N` nodes.


  n = nrow(TS)
  W = matrix(0, n, n)

  for(i in 1 : n){
    xi = split_data(TS[i, ], lag)$inputs
    yi = split_data(TS[i, ], lag)$targets

    for (j in 1:n) {
      xj = split_data(TS[j, ], lag)$inputs
      yj = split_data(TS[j, ], lag)$targets
      xij = cbind(xi, xj)
      reg1 = lm(as.vector(yi) ~  xi)
      reg2 = lm(as.vector(yi) ~  xij)
      err1 = yi - predict(reg1, as.data.frame(xi))
      err2 = suppressWarnings(yi - predict(reg2, as.data.frame(xij)))

      std_i = pracma::std(err1, 1)
      std_ij = pracma::std(err2, 1)

      if(std_i == 0){
        W[j, i] = -99999999
      } else if(std_ij == 0){
        W[j, i] = 99999999
      } else {
        W[j, i] = log(std_i) - log(std_ij)
      }
    }
  }

  # threshold the network
  W_thresh = threshold(W, threshold_type, ...)

  # construct the network
  G = create_graph(W_thresh)

  structure(
    list(
      weights_matrix = W,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "GrangerCausality"
  )
}


split_data <- function(TS, lag){
        # From a single node time series, return a training dataset with
        # corresponding targets.
        #
        # Parameters
        # ----------
        #
        # TS
        #     Original code state it shoudl be an array consisting of :math:`L`
        # observations from :math:`N` sensors.
        #     However, it should be an 1-D Matrix
        #
        # lag (int)
        #     Time lag to consider.
        #
        # Returns
        # -------
        #
        # inputs
        #     Training data for the inputs.
        #
        # targets
        #     Training data for the targets.

  T = length(TS)
  inputs = matrix(0, (T - lag - 1), lag)
  targets = matrix(0, (T - lag - 1))

  for (t in 1:((T - lag - 1))) {
    inputs[t,] = TS[lag + t - 1]
    targets[t] = TS[t + lag]
  }

  mylist <- list("inputs" = inputs, "targets" = t(targets))
  return(mylist)
}



#test
# TS = matrix(c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9), nrow=3, ncol=7, byrow = TRUE)
# Ts = TS[2,]
# lag = 1
# # print(split_data(Ts, lag)$inputs)
# # print(split_data(Ts, lag)$targets)
# s <- granger_causality_fit(TS)
# s
# plot(s$graph)











