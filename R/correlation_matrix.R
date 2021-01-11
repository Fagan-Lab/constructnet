#correlation_matrix.R
#status Unfinished
#Reconstruction of graphs using the correlation matrix.
#author Zhaoyi Zhuang

fit <- function(TS, num_eigs=Null, threshold_type='range', ...) {

        #If ``num_eigs`` is `Null`, perform the reconstruction using the
        #unregularized correlation matrix. Otherwise, construct a regularized
        #precision matrix using ``num_eigs`` eigenvectors and eigenvalues of the
        #correlation matrix.

        #Parameters
        #----------
        #TS matrix
        #    Matrix consisting of :`L` observations from :`N` sensors

        #num_eigs (int)
        #    The number of eigenvalues to use. (This corresponds to the
        #    amount of regularization.) The number of eigenvalues used must
        #    be less than :`N`.

        #threshold_type (str)
        #    Which thresholding function to use on the matrix of
        #    weights.

  # return list
  results <- list()

  # get the correlation matrix
  co = cor(t(TS), method="pearson")

  if (!is.null(num_eigs)){
    N <- ncol(TS)
    if (num_eigs >  N) {
      stop("The number of eigenvalues used must be less than the number of sensors.")
    }

    # get eigenvalues and eigenvectors of the correlation matrix
    ev <- eigen(t(co))
    values <- ev$values
    vectors <- ev$vectors

    # construct the precision matrix and store it
    P = (vectors[,c(1:num_eigs)]) %*%
        (c (1 / matrix(values[1:num_eigs], ncol = 1, byrow = TRUE)) * t(vectors[,c(1:num_eigs)]))
    P = P / (
        matrix(sqrt(diag(P)), ncol = 1, byrow = TRUE) %*% matrix(sqrt(diag(P)), nrow = 1, byrow = TRUE)
    )
    mat <- P
  } else {
    mat <- co
  }

  #waiting for threshold.R
  #A = threshold(mat, threshold_type, ...)
  A <- matrix(c(2, 4, 3, 1, 5, 7), nrow=2, ncol=3, byrow = TRUE)


  structure(
    list(
      # store the appropriate source matrix
      weights_matrix = mat,
      # threshold the correlation matrix
      thresholded_matrix = A,
      # construct the network
      # waiting for graph.R
      graph = A
      #graph = create_graph(A)
    ),
    class = "correlationMatrix"
  )

}


# test
#mydata = matrix(c(2, 4, 3, 1, 5, 7, 3, 4, 7), nrow=3, ncol=3, byrow = TRUE)
#fit(mydata, 2)

