# maximum_likelihood_estimation.R
# ---------------------
# status: Finished draft and some simple tests
# Reconstruction of graphs using maximum likelihood estimation
# author: Brennan Klein
# converted by: Zhaoyi Zhuang

devtools::load_all(".")
# source(here::here('constructnet', 'R', 'threshold.R'))
# source(here::here('constructnet', 'R', 'graph.R'))

maximum_likelihood_estimation_fit <- function(TS, rate=1.0, stop_criterion=T, threshold_type='degree', ...){
    # Infer inter-node coupling weights using maximum likelihood estimation
    # methods.
    # The results dictionary also stores the weight matrix as
    # `'weights_matrix'` and the thresholded version of the weight matrix
    # as `'thresholded_matrix'`.
    #
    # Parameters
    # ----------
    #
    # TS
    #   Matrix consisting of :math:`L` observations from :math:`N` sensors.
    #
    # rate (float)
    #   rate term in maximum likelihood
    #
    # stop_criterion (bool)
    #    if True, prevent overly-long runtimes
    #
    # threshold_type (str)
    #    Which thresholding function to use on the matrix of
    # weights.
    #
    # Returns
    # -------
    #   G
    #      a reconstructed graph.
  N = nrow(TS)
  L = ncol(TS)
  rate = rate/L
  W = matrix(0, N, N)




}



#test






