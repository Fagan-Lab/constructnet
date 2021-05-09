#' Reconstruction of graphs using the partial correlation influence.
#'
#'
#' @param TS Input matrix, X * L matrix consisting of L observations from N seneors.
#' @param index An index variable or set of index variables, which are assumed
#' to be confounders of all other variables.
#' @param threshold_type Which thresholding function to use on the matrix of weights
#' @param ... arguments
#'
#' @return A reconstructed graph with N nodes
#' @export
partial_correlation_influence <- function(TS, index=None, threshold_type='range', ...) {
  data <- t(TS)
  N = nrow(data)

  mask = matrix(1, ncol = N)
  if (!is.null(index)) {
    mask[index] <- FALSE
  }
  j <- 1;
  mask_index <- vector()
  for (i in 1:N) {
    if (mask[i] == 1) {
      mask_index[j] = i
      j = j + 1
    }
  }
  
  p_corr = matrix(NA, N, N)
  p_corr[mask_index, mask_index] = partial_corr(data[, mask_index], data[, !mask])
  p_corr = t(p_corr)
  
  # For every non-index variable Z, compute partial correlation influence
  # between other variables when Z is also held constant
  p_corr_inf <- list()
  for(i in 1:N) {
    p_corr_inf[[i]] <- matrix(NA, N, N)
  }
  
  for (z in seq(N)[mask_index]) {
    m_new = mask
    m_new[z] = FALSE
    
    j <- 1;
    m_index <- vector()
    for (i in 1:N) {
      if (m_new[i] == 1) {
        m_index[j] = i 
        j = j + 1
      }
    }
    
    diff = p_corr[m_index, m_index]
    diff = diff - partial_corr(data[, m_index], data[, !m_new])
    
    # p_corr_inf[np.ix_(m_new, m_new, [z])] = diff[:, :, np.newaxis]
    
    diff_values = c(t(diff))
    i <- 1
    for(x in m_index) {
      for(c in m_index) {
        p_corr_inf[[x]][c,z] = diff_values[i]
        i = (i+1)
        if(i > length(diff_values)) {
          i = i%%length(diff_values)
        }
      }
    }
    
    # np.fill_diagonal(p_corr_inf[:, :, z], np.nan)
    for(i in 1:N) {
      p_corr_inf[[i]][i,N] = NA
    }
    
    p_corr_inf[[z]][, z] = 0
  }
  
  influence <- matrix(0, N, N)
  # influence[mask, mask] = np.nanmean(p_corr_inf[mask, mask], axis=1)
  for(i in 1:N) {
    if (mask[i] == 1) {
      influence[i,i] = mean(p_corr_inf[[i]][i,], na.rm=TRUE)
    }
  }
    
  influence[!mask, ] = Inf
  influence[, !mask] = 0
  
  W_thresh = threshold(influene, threshold_type, ...)
  G = create_graph(W_thresh)
  
  structure(
    list(
      weights_matrix = influence,
      thresholded_matrix = W_thresh,
      graph = G
    ),
    class = "PartialCorrelationInfluence"
  )
} 

partial_corr <- function(vars, idx_vars) {
  if(length(idx_vars) == 0) {
    return(stats::cor(t(vars), method = "pearson"))
  } else {
    coef = solve(vars, idx_vars)[1]
    resid = vars - (idx_vars %*% coef)
    return(stats::cor(t(resid), method = "pearson"))
  }
}