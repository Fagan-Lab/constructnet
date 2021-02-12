
mask_function <- function(mat, cutoffs) {
  # setting values not within a list of ranges to zero and others to one.

  # Parameters
  #----------
  # mat
  #   matrix

  # cutoffs (list of tuples)
  #    When thresholding, include only edges whose correlations fall
  #    within a given range or set of ranges. The lower value must come
  #    first in each tuple.
  # Returns
  #-------
  #  mask
  #     matrix
  mask <- mat
  for (row in 1:nrow(mat)) {
    for (col in 1:ncol(mat)) {
      for (cutoff in cutoffs) {
        mask[row, col] <- (mat[row, col] >= cutoff[1] && mat[row, col] <= cutoff[2])
      }
    }
  }
  result <- mask
}

threshold_in_range <- function(mat, ...) {
  # Threshold by setting values not within a list of ranges to zero.

  # Parameters
  #----------
  # mat
  #   matrix

  # cutoffs (list of tuples)
  #    When thresholding, include only edges whose correlations fall
  #    within a given range or set of ranges. The lower value must come
  #    first in each tuple.
  # Returns
  #-------
  #  thresholded_mat
  #     the thresholded matrix


  kwargs <- list(...)
  if ("cutoffs" %in% names(kwargs)) {
    cutoffs <- kwargs[["cutoffs"]]
  } else {
    warning("Setting 'cutoffs' argument is strongly encouraged. Using cutoff range of (-1, 1).")
    cutoffs <- list(c(-1, 1))
  }

  mask <- mask_function(mat, cutoffs)
  thresholded_mat <- mat * mask

  if (!is.null(kwargs[["binary"]])) {
    if (kwargs[["binary"]] || FALSE) {
      thresholded_mat <- mask
    }
  }

  if (!is.null(kwargs[["remove_self_loops"]])) {
    if (kwargs[["remove_self_loops"]] && TRUE) {
      diag(thresholded_mat) <- 0
    }
  } else {
    diag(thresholded_mat) <- 0
  }

  result <- thresholded_mat
}


threshold_on_quantile <- function(mat, ...) {
  # Threshold by setting values below a given quantile to zero.

  # Parameters
  #----------

  # mat
  #   matrix

  # quantile (float)
  #    The threshold above which to keep an element of the array, e.g.,
  #    set to zero elements below the 90th quantile of the array.

  # Returns
  #-------
  # thresholded_mat
  #    the thresholded matrix

  kwargs <- list(...)
  if ("quantile" %in% names(kwargs)) {
    quantile <- kwargs[["quantile"]]
  } else {
    warning("Setting 'quantile' argument is strongly recommended. Using target quantile of 0.9 for thresholding.")
    quantile <- 0.9
  }

  if (!is.null(kwargs[["remove_self_loops"]])) {
    if (kwargs[["remove_self_loops"]] && TRUE) {
      diag(mat) <- 0
    }
  } else {
    diag(mat) <- 0
  }

  if (quantile != 0) {
    thresholded_mat <- mat * (mat > quantile(mat, probs = quantile))
  } else {
    thresholded_mat <- mat
  }

  if (!is.null(kwargs[["binary"]])) {
    if (kwargs[["binary"]] || FALSE) {
      thresholded_mat <- abs(sign(thresholded_mat))
    }
  }

  result <- thresholded_mat
}


threshold_on_degree <- function(mat, ...) {
  # Threshold by setting values below a given quantile to zero.

  # Parameters
  #----------

  # mat
  #  matrix

  # avg_k (float)
  #    The average degree to target when thresholding the matrix.

  # Returns
  #-------
  # thresholded_mat
  #    the thresholded matrix
  kwargs <- list(...)
  if ("avg_k" %in% names(kwargs)) {
    avg_k <- kwargs[["avg_k"]]
  } else {
    warning("Setting 'avg_k' argument is strongly encouraged. Using average degree of 1 for thresholding.")
    avg_k <- 1
  }

  n <- ncol(mat)
  A <- matrix(1, n, n)

  if (!is.null(kwargs[["remove_self_loops"]])) {
    if (kwargs[["remove_self_loops"]] && TRUE) {
      diag(A) <- 0
      diag(mat) <- 0
    }
  } else {
    diag(A) <- 0
    diag(mat) <- 0
  }

  if (mean(rowSums(A)) <= avg_k) {
    thresholded_mat <- mat
  } else {
    for (m in sort(as.vector(t(mat)))) {
      A[mat == m] <- 0
      if (mean(rowSums(A)) <= avg_k) {
        break
      }
    }
    thresholded_mat <- mat * (mat > m)
  }

  if (!is.null(kwargs[["binary"]])) {
    if (kwargs[["binary"]] || FALSE) {
      thresholded_mat <- abs(sign(thresholded_mat))
    }
  }

  result <- thresholded_mat
}


#' Utilities for thresholding matrices based on different criteria
#' @param mat input matrix
#'
#' @param rule A string indicating which thresholding function to invoke.
#' @param ... Arguments
#'
#' @export
threshold <- function(mat, rule, ...) {
  # A flexible interface to other thresholding functions.

  # Returns
  #-------
  # thresholded_mat
  #    the thresholded matrix

  kwargs <- list(...)
  result <- tryCatch(
    {
      if (rule == "degree") {
        threshold_on_degree(mat, ...)
      }
      else if (rule == "range") {
        threshold_in_range(mat, ...)
      }
      else if (rule == "quantile") {
        threshold_on_quantile(mat, ...)
      }
      else if (rule == "custom") {
        kwargs[["custom_thresholder"]](mat)
      }
    },
    error = function(e) {
      stop("missing threshold parameter")
    }
  )
}
