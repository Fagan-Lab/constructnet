# marchenko_pastur.R
# --------------
# status: Finished draft and some simple tests
# Graph reconstruction algorithm.
# author: Matteo Chinazzi
# converted by: Zhaoyi Zhuang

library(igraph)
source(here::here('constructnet', 'R', 'threshold.R'))
source(here::here('constructnet', 'R', 'graph.R'))

marchenko_pastur_fit <- function(TS, remove_largest=F, metric_distance=F, tol = 1e-15, threshold_type='range', ...){
  # Create a correlation-based graph using Marchenko-Pastur law to remove noise.
  # Parameters
  # ----------
  #   
  # TS 
  #   N times L matrix consisting of `L` observations
  #   from `N` sensors.
  # 
  # remove_largest (bool)
  #
  # tol
  #  avoid minor difference between python and r in case of calculating eigen values
  # 
  # Returns
  # -------
  #   
  #   G 
  #     A reconstructed graph.
  N = nrow(TS)
  L = ncol(TS)
  if(N > L){
    stop('"L must be greater or equal than N."')
  }
  
  Q = L/N
  C = cor(t(TS))
  ev <- eigen(t(C))
  w <- ev$values
  v <- ev$vectors
  w_min = 1 + 1 / Q - 2 * sqrt(1 / Q)
  w_max = 1 + 1 / Q +  2 * sqrt(1 / Q)
  
  selected = (w < w_min+tol) | (w > w_max-tol)
  if(sum(selected) == 0){
    G = make_empty_graph(n = N)
    result <- G
    return(result)
  }

  if(remove_largest){
    selected[1] <- F
  }

  w_signal = as.matrix(w[selected])
  v_signal = as.matrix(v[, selected])
  C_signal = v_signal %*% (diag(w_signal)) %*% t(v_signal)
  
  if(metric_distance){
    C_signal = sqrt(2 * ( 1 - C_signal))
  }
  
  t = threshold(C_signal, threshold_type, ...)
  
  G = create_graph(t)
  
  structure(
    list(
      weights_matrix = C_signal,
      thresholded_matrix = t,
      graph = G
    ),
    class = "GraphicalLasso"
  )
  
  print(C_signal)
  print(t)
  result <- G
  return(result)
}







#test
# TS = matrix(c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9), nrow=3, ncol=7, byrow = TRUE)
# m <- marchenko_pastur_fit(TS, T)
# print(m)
# plot(m)
# 
# TS = matrix(c(6, 3, 1, 5, 3, 0, 5, 1, 1, 5, 6, 2, 5, 1, 2, 2), nrow = 4, ncol = 4)
# m <- marchenko_pastur_fit(TS, T, T)
# plot(m)






