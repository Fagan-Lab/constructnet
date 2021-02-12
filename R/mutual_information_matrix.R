#' Calculates the mutual information between the probability distributions
#' of the (binned) values of the time series of pairs of nodes.
#' First, the mutual information is computed between each pair of
#' vertices.  Then, a thresholding condition is applied to obtain
#' edges.
#' The results dictionary also stores the weight matrix as
#' `'weights_matrix'` and the thresholded version of the weight matrix
#' as `'thresholded_matrix'`.
#'
#' @param TS Matrix consisting of :math:`L` observations from :math:`N` sensors.
#' @param nbins number of bins for the pre-processing step (to yield a discrete probability distribution)
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param ... arguments
#'
#' @return A reconstructed graph with :math:`N` nodes.
#' @export

mutual_information_matrix_fit <- function(TS, nbins = 10, threshold_type = "degree", ...) {
  N <- nrow(TS)
  rang <- c(min(TS), max(TS))

  IndivP <- find_individual_probability_distribution(TS, rang, nbins)
  ProduP <- find_product_probability_distribution(IndivP, N)
  JointP <- find_joint_probability_distribution(TS, rang, nbins)

  I <- mutual_info_all_pairs(JointP, ProduP, N)

  A <- threshold(I, threshold_type, ...)
  G <- create_graph(A)

  structure(
    list(
      weights_matrix = I,
      thresholded_matrix = A,
      graph = G
    ),
    class = "MutualInformationMatrix"
  )
}

hist2d <- function(x,
                   y = NULL,
                   nbins = 200,
                   custom,
                   same.scale = FALSE,
                   na.rm = TRUE,
                   show = TRUE,
                   col = c("black", grDevices::heat.colors(12)),
                   FUN = base::length,
                   xlab,
                   ylab,
                   ...) {

  # custom hist2d.
  # add designed nbins and change include of right to exclude.
  if (is.null(y)) {
    if (ncol(x) != 2) stop("If y is ommitted, x must be a 2 column matirx")
    y <- x[, 2]
    x <- x[, 1]
  }

  if (length(nbins) == 1) {
    nbins <- rep(nbins, 2)
  }

  nas <- is.na(x) | is.na(y)

  if (na.rm) {
    x <- x[!nas]
    y <- y[!nas]
  }
  else {
    stop("missinig values not permitted if na.rm=FALSE")
  }

  if (same.scale) {
    x.cuts <- seq(from = min(x, y), to = max(x, y), length = nbins[1] + 1)
    y.cuts <- seq(from = min(x, y), to = max(x, y), length = nbins[2] + 1)
  }
  else {
    x.cuts <- custom
    y.cuts <- custom
  }

  index.x <- cut(x, x.cuts, right = F, include.lowest = T)
  index.y <- cut(y, y.cuts, right = F, include.lowest = T)

  ## tapply is faster than old for() loop, and allows
  ## use of any user-specified summary function
  m <- tapply(x, list(index.x, index.y), FUN)

  ## If we're using length, set empty cells to 0 instead of NA
  if (identical(FUN, base::length)) {
    m[is.na(m)] <- 0
  }

  if (missing(xlab)) xlab <- deparse(substitute(xlab))
  if (missing(ylab)) ylab <- deparse(substitute(ylab))

  if (show) {
    graphics::image(x.cuts, y.cuts, m, col = col, xlab = xlab, ylab = ylab, ...)
  }

  midpoints <- function(x) (x[-1] + x[-length(x)]) / 2

  retval <- list()
  retval$counts <- m
  retval$x.breaks <- x.cuts
  retval$y.breaks <- y.cuts
  retval$x <- midpoints(x.cuts)
  retval$y <- midpoints(y.cuts)
  retval$nobs <- length(x)
  retval$call <- match.call()
  class(retval) <- "hist2d"
  retval
}

find_individual_probability_distribution <- function(TS, rang, nbins) {
  N <- nrow(TS)
  L <- ncol(TS)
  IndivP <- Dict::dict()

  x <- seq(from = rang[1], to = rang[2], length = nbins + 1)

  for (j in 1:N) {
    s <- graphics::hist(TS[j, ], breaks = x, right = F, include.lowest = T)$counts
    IndivP[toString(j)] <- s / L
  }

  IndivP
}

find_product_probability_distribution <- function(IndivP, N) {
  ProduP <- Dict::dict()

  for (l in 1:N) {
    for (j in 1:l) {
      if (l != j) {
        ProduP[toString(c(j, l))] <- IndivP[j] %o% IndivP[l]
      }
    }
  }
  ProduP
}

find_joint_probability_distribution <- function(TS, rang, nbins) {
  N <- nrow(TS)
  L <- ncol(TS)
  JointP <- Dict::dict()

  x <- seq(from = rang[1], to = rang[2], length = nbins + 1)


  for (l in 1:N) {
    for (j in 1:l) {
      if (l != j) {
        s <- hist2d(TS[j, ], TS[l, ], nbins = nbins, custom = x)$counts
        JointP[toString(c(j, l))] <- s / L
      }
    }
  }
  JointP
}


mutual_info_node_pair <- function(JointP_jl, ProduP_jl) {
  I_jl <- 0

  z <- data.frame(x = as.vector(t(JointP_jl)), y = as.vector(t(ProduP_jl)))

  for (i in 1:nrow(z)) {
    row <- z[i, ]
    if (row$x > 0 && row$y > 0) {
      I_jl <- I_jl + row$x * log(row$x / row$y)
    }
  }

  I_jl
}

mutual_info_all_pairs <- function(JointP, ProduP, N) {
  I <- matrix(0, N, N)

  for (l in 1:N) {
    for (j in 1:l) {
      if (l != j) {
        JointP_jl <- JointP[toString(c(j, l))]
        ProduP_jl <- ProduP[toString(c(j, l))]

        I[j, l] <- mutual_info_node_pair(JointP_jl, ProduP_jl)
        I[l, j] <- I[j, l]
      }
    }
  }
  I
}

threshold_from_degree <- function(deg, M) {
  N <- nrow(M)
  A <- matrix(1, N, N)
  for (tau in sort(as.vector(M))) {
    idx <- which(M == tau)
    for (i in idx) {
      A[i] <- 0
    }
    if (mean(rowSums(A)) < deg) {
      break
    }
  }
  return(tau)
}
