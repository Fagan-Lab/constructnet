threshold_from_degree <- function(deg, M){

  # Compute the required threshold (tau) in order to yield a reconstructed graph of mean degree deg.
  # Parameters
  # ----------
  #   deg (int): Target degree for which the appropriate threshold will be computed
  #   M : Pre-thresholded NxN array
  # Returns
  # ------
  #   tau (float): Required threshold for A=np.array(I<tau,dtype=int) to have an average of deg ones per row/column

N = nrow(M)
A = matrix(1, N, N)
for(tau in sort(as.vector(M))) {
  idx = which(M == tau)
  for(i in idx){
    A[i] <- 0
  }
  if(mean(rowSums(A)) < deg){
    break
  }
}

return(tau)

}



#test
TS = matrix(c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9), nrow=3, ncol=7, byrow = TRUE)
M = matrix(c(5, 1, 16, 2, -5, -1, 0, 2, 14, 2, 61, 6, 13 , 14 , 8, 9), nrow=4, ncol=4, byrow = TRUE)
deg = 1
threshold_from_degree(deg, M)

