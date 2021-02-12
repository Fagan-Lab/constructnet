
#' Reconstruction of graphs using the correlation matrix.

#' @param TS Matrix consisting of L observations from N sensors
#'
#' @param num_eigs The number of eigenvalues to use. (This corresponds to the
#'                 amount of regularization.) The number of eigenvalues used must
#'                 be less than N.
#'
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#'
#' @param ... ARGUMENTS
#'
#' @export
correlation_matrix_fit <- function(TS, num_eigs = NULL, threshold_type = "range", ...) {

  # If ``num_eigs`` is `Null`, perform the reconstruction using the
  # unregularized correlation matrix. Otherwise, construct a regularized
  # precision matrix using ``num_eigs`` eigenvectors and eigenvalues of the
  # correlation matrix.

  # Returns
  #-------
  #  G
  #   a reconstructed graph.
  #  list()
  #   structure


  # get the correlation matrix
  co <- cor(t(TS), method = "pearson")

  if (!is.null(num_eigs)) {
    N <- ncol(TS)
    if (num_eigs > N) {
      stop("The number of eigenvalues used must be less than the number of sensors.")
    }

    # get eigenvalues and eigenvectors of the correlation matrix
    ev <- eigen(t(co))
    values <- ev$values
    vectors <- ev$vectors

    # construct the precision matrix and store it
    P <- (vectors[, c(1:num_eigs)]) %*%
      (c(1 / matrix(values[1:num_eigs], ncol = 1, byrow = TRUE)) * t(vectors[, c(1:num_eigs)]))
    P <- P / (
      matrix(sqrt(diag(P)), ncol = 1, byrow = TRUE) %*% matrix(sqrt(diag(P)), nrow = 1, byrow = TRUE)
    )
    mat <- P
  } else {
    mat <- co
  }

  A <- threshold(mat, threshold_type, ...)

  structure(
    list(
      # store the appropriate source matrix
      weights_matrix = mat,
      # threshold the correlation matrix
      thresholded_matrix = A,
      # construct the network
      graph = create_graph(A)
    ),
    class = "correlationMatrix"
  )
}
