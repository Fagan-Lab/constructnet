# free_energy_minimization.R
# ---------------------------
# status: Finished draft and some simple tests
# Reconstruction of graphs by minimizing a free energy of your data
# author: Brennan Klein
# converted by: Zhaoyi Zhuang

source(here::here('constructnet', 'R', 'threshold.R'))
source(here::here('constructnet', 'R', 'graph.R'))

free_energy_minimization_fit <- function(TS, threshold_type='degree', ...){
        # Infer inter-node coupling weights by minimizing a free energy over the
        # data structure.
        # 
        # Parameters
        # ----------
        # 
        # TS 
        #     Matrix consisting of :`L` observations from :`N`
        #     sensors.
        # 
        # threshold_type (str)
        #     Which thresholding function to use on the matrix of
        #     weights.
        # 
        # Returns
        # -------
        # 
        # G 
        #     a reconstructed graph.
  N = nrow(TS)
  L = ncol(TS)
  m <- TS[, -L]
  m <- rowMeans(m) # model average
  ds <- t(TS[, -L])
  ds <- apply(ds, 1, function(x) x - m)
  ds <- t(ds)      # discrepancy
  t1 = L -1        # time limit
  
  # covariance of the discrepeancy
  c = cov.wt(ds, method = 'ML')$cov
  c_inv = solve(c)
  dst = t(ds)       # discrepancy at time t
  W = matrix(1, N, N)
  nloop = 3     # failsafe
  
  for(i0 in 1:N){
    TS1 = TS[i0, -1]
    h = TS1
    cost = matrix(100, 1, nloop)
    for(iloop in 1:nloop){
      h_av = mean(h)                        # average local field
      hs_av = t(dst %*% (h - h_av) / t1)  # deltaE_i delta\sigma_k 
      w = (hs_av %*% c_inv)                 # expectation under model
      h = as.vector(t(t(TS[, -L]) %*% w[,]))  # estimate of local field
      TS_model = tanh(h)                     # under kinetic Ising model
      # discrepancy cost
      cost[iloop] = mean((TS1 - TS_model)**2 )
      
      if(iloop != 1 && (cost[iloop] >= cost[iloop - 1])) {
        break
      }
      
      #The original Python code does not work 
      #
      # h *= np.divide(
      #   TS1, TS_model, out=np.ones_like(TS1), where=TS_model != 0
      # )
      #
      #It works when I delete "out=np.ones_like(TS1)"
      h = h *  (TS1 / TS_model)
    }
    W[i0, ] = w[,]
  }
  
  # threshold the network
  W_thresh = threshold(W, threshold_type, ...)
  
  
  G = create_graph(W_thresh)
  # construct the network
  structure(
    list(
      weights_matrix = W,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "FreeEnergyMinimization"
  )

  # in the original Python code, function threshold() has side effect 
  # which change the value of W to W_thresh. In R, W will not be changed
  # even after threshold() is applied based on it. 
  print(W) 
  print(W_thresh)
  results <- G
}

#TEST
# TS = matrix(c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9), nrow=3, ncol=7, byrow = TRUE)
# f <- free_energy_minimization_fit(TS)
# print(f)
# plot(f)
