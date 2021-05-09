#' Reconstruction of graphs using the partial correlation matrix.
#'
#'
#' @param TS Input matrix, X * L matrix consisting of L observations from N seneors.
#' @param index An index variable or set of index variables, which are assumed
#' to be confounders of all other variables.
#' @param drop_index If True, drop the index variables after calculating the partial
#' correlations.
#' @param of_residuals If True, after calculating the partial correlations (presumably
#' using a dropped index variable), recalculate the partial correlations between 
#' each variable, holding constant all other variables.
#' @param threshold_type Which thresholding function to use on the matrix of weights
#' @param ... arguments
#'
#' @return A reconstructed graph with N nodes
#' @export
partial_correlation_matrix <- function(TS, index=NULL, drop_index=TRUE, 
                                       of_residuals=FALSE, threshold_type='range', ...) {
  
  p_cor = partial_corr(TS, index=index)
  
  if(drop_index && !is.null(index)) {
    p_cor = p_cor[-index,]
    p_cor = p_corr[, -index]
  }
  
  if(of_residuals) {
    p_cor = partial_corr(p_cor, index=NULL)
  }
  
  W_thresh = threshold(p_cor, threshold_type, ...)
  G = create_graph(W_thresh)
  
  structure(
    list(
      weights_matrix = p_cor,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "PartialCorrelationMatrix"
  )
}

partial_corr <- function(C, index=NULL) {
  C = t(C)
  p <- nrow(C)
  P_corr <- matrix(0.0, N, N)
  
  for(i in 1:p) {
    P_corr[i,i] = 1
    if(i < p){
      for(j in (i+1):p){
        if(!is.null(index)){
          idx <- matrix(TRUE, ncol= p)
          idx[i] = FALSE
          idx[j] = FALSE
        }else if(typeof(index)=="double" || (is(index, "array") && mode(tst) == "numeric") ){
          idx <- matrix(FALSE, ncol= p)
          idx[index] = TRUE
        }else{
          stop("Index must be an integer, an array of " "integers, or None.")
        }
        
        beta_i = solve(C[, idx], C[, j])[1]
        beta_j = solve(C[, idx], C[, i])[1]
        
        res_j = C[, j] - (C[, idx] %*% beta_i)
        res_i = C[, i] - (C[, idx] %*% beta_j)
        
        corr = stats::cor.test(res_i, res_j)[1]
        P_corr[i, j] = corr
        P_corr[j, i] = corr
      }
    }
  }
  
  P_corr
}