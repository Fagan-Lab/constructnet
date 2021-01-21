# graphical_lasso.R
# --------------
# status: Finished draft and some simple tests
# Graph reconstruction algorithm.
# author: Charles Murphy
# converted by: Zhaoyi Zhuang

# source(here::here('constructnet', 'R', 'threshold.R'))
# source(here::here('constructnet', 'R', 'graph.R'))

graphical_lasso_fit <- function(TS, alpha = 0.01, max_iter = 100, tol = 0.0001, threshold_type = 'degree', ...){
        # Performs a graphical lasso.
        #
        # Parameters
        # ----------
        #
        # TS
        #     Matrix consisting of :`L` observations from :`N`
        #     sensors.
        #
        # alpha
        #     Coefficient of penalization, higher values means more
        #     sparseness
        #
        # max_iter
        #     Maximum number of iterations.
        #
        # tol
        #     Stop the algorithm when the duality gap is below a certain
        #     threshold.
        #
        # threshold_type
        #     Which thresholding function to use on the matrix of
        #     weights.
        #
        # Returns
        # -------
        #
        # G
        #     A reconstructed graph with :`N` nodes.

  emp_cov = cov(t(TS))
  prec = glasso::glasso(emp_cov, alpha, maxit = max_iter, thr = tol)$wi
  cov = glasso::glasso(emp_cov, alpha, maxit = max_iter, thr = tol)$w

  W_thresh = threshold(cov, threshold_type, ...)

  G = create_graph(W_thresh)

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



#test
# TS = matrix(c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9), nrow=3, ncol=7, byrow = TRUE)
# x <- graphical_lasso_fit(TS)
# x
# plot(x$graph)





