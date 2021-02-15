joint_entropy <-  function(data){
   # Joint entropy of all variables in the data.
   # Parameters
   # ----------
   #   data
   # matrix of data with variables as columns and observations as rows.
   #
   #  Returns
   #  -------
   #  float
   #      Joint entropy of the variables of interests.

  count = Dict::dict()
  N = nrow(data)
  L = ncol(data)

  c = matrix(0, N, 1)
  if(L == 1) {
    data <- cbind(data, c)
    data <- cbind(data, c)
  }

  if(L == 2) {
    data <- cbind(data, c)
  }

  for (i in 1:N) {
    count[toString(c(data[i, 1], data[i, 2], data[i, 3]))] = 0
  }
  for (i in 1:N) {
    count[toString(c(data[i, 1], data[i, 2], data[i, 3]))] = count[toString(c(data[i, 1], data[i, 2], data[i, 3]))] + 1
  }

  h = unlist(count$values)
  DescTools::Entropy(h , base = 2)
}

#' Conditional entropy of variables in the data conditioned on
#' a given set of variables.
#' @param data matrix of data with variables of interests as columns and
#' observations as rows.
#' @param given matrix of data with the conditioned variables as columns and
#' observations as rows.
#'
#' @export
#'
conditional_entropy <- function(data, given){

   #
   #  Returns
   #  -------
   #  float
   #      Conditional entrpoy of the variables X_i of interest
   #      conditioned on variables Y_j.
  joint = cbind(data, given)
  entrp = joint_entropy(joint) - joint_entropy(given)
  entrp
}

#' An entry in the returned array is the index of the bin of the
#' linearly-binned raw continuous data.
#' @param raw matrix of raw continuous data.
#' @param n_bins A universal number of bins for all the variables.
#'
#' @export
#'
categorized_data <- function(raw, n_bins){
  # Returns
  # -------
  # Matrix of bin indices after categorizing the raw data.
  N = nrow(raw)
  L = ncol(raw)
  bins = linear_bins(raw, n_bins)
  data = matrix(1, N, L)
  for(i in 1:N){
    for (j in 1:L){
      x = t(as.matrix(bins[1:n_bins+1, j]>=raw[i, j]) %*% 1)
      data[i,j] <- ramify::argmax(x)
    }
  }
  data
}


linear_bins <- function(raw, n_bins){
  # Separators of linear bins for each variable in the raw data.
  #
  # Parameters
  # ----------
  #   raw
  # matrix of raw continuous data.
  #
  # n_bins (int)
  # A universal number of bins for all the variables.
  #
  # Returns
  # -------
  # matrix where a column is the separators of bins for a variable.
  minm = apply(raw, 2, min)
  maxm = apply(raw, 2, max)
  bins <- vector()
  for (i in 1:length(minm)) {
    bins <- append(bins, pracma::linspace(minm[i], maxm[i], n = n_bins+1))
  }
  bins <- matrix(bins, nrow = n_bins + 1, ncol = length(minm), byrow = F)
  bins
}
