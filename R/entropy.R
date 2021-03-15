joint_entropy <- function(data) {
  count <- Dict::dict()
  N <- nrow(data)
  L <- ncol(data)
  
  c <- matrix(0, N, 1)
  if (L == 1) {
    data <- cbind(data, c)
  }
  
  for (i in 1:N) {
    count[toString(c(data[i, 1], data[i, 2], data[i, 3]))] = 0
  }
  for (i in 1:N) {
    count[toString(c(data[i, 1], data[i, 2], data[i, 3]))] = count[toString(c(data[i, 1], data[i, 2], data[i, 3]))] + 1
  }
  
  h <- unlist(count$values)
  DescTools::Entropy(h, base = 2)
}

#' Conditional entropy of variables in the data conditioned on
#' a given set of variables.
#'
#' @param data matrix of data with variables of interests as columns and
#' observations as rows.
#' @param given matrix of data with the conditioned variables as columns and
#' observations as rows.
#'
#' @export
#'
conditional_entropy <- function(data, given) {
  joint <- cbind(data, given)
  entrp <- joint_entropy(joint) - joint_entropy(given)
  entrp
}

#' An entry in the returned array is the index of the bin of the
#' linearly-binned raw continuous data.
#' @param raw matrix of raw continuous data.
#' @param n_bins A universal number of bins for all the variables.
#'
#' @export
#'
categorized_data <- function(raw, n_bins) {
  N <- nrow(raw)
  L <- ncol(raw)
  bins <- linear_bins(raw, n_bins)
  data <- matrix(1, N, L)
  for (i in 1:N) {
    for (j in 1:L) {
      x <- t(as.matrix(bins[1:n_bins + 1, j] >= raw[i, j]) %*% 1)
      data[i, j] <- ramify::argmax(x)
    }
  }
  
  data
}


linear_bins <- function(raw, n_bins) {
  minm <- apply(raw, 2, min)
  maxm <- apply(raw, 2, max)
  bins <- vector()
  
  for (i in 1:length(minm)) {
    bins <- append(bins, pracma::linspace(minm[i], maxm[i], n = n_bins + 1))
  }
  
  matrix(bins, nrow = n_bins + 1, ncol = length(minm), byrow = F)
}