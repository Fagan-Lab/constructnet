
#' Infer inter-node coupling weights using a Thouless-Anderson-Palmer mean
#' field approximation.
#'
#' @param TS Matrix consisting of L observations from N sensors.
#' @param threshold_type Which thresholding function to use on the matrix of weights.
#' @param ... arguments
#'
#' @return a reconstructed graph.
#' @export
thouless_anderson_palmer_fit <- function(TS, threshold_type='range',...){
    N = nrow(TS)
    L = ncol(TS)
    m = rowMeans(TS)

    A = 1 - m ** 2
    A_inv = diag(1/A)
    A = diag(A)
    ds = sweep(t(TS) , 2, m, '-')
    C = cov.wt(ds, method = 'ML')$cov
    C_inv = solve(C)

    s1 = TS[, -1]
    ds1 = sweep(t(s1) , 2, rowMeans(s1), '-')
    D = cross_cov(ds1, ds[-(nrow(ds)),])

    B = D %*% C_inv
    W_NMF = A_inv %*% B

    step = 0.001
    nloop = as.integer(0.33/step) + 2

    W2_NMF = W_NMF **2
    temp = matrix(0, 1, N)
    F1 = matrix(0, 1, N)

    for (i in 1:N) {
      temp[i] = (1 - m[i]**2) * sum(W2_NMF[i,] * (1-m**2))
      y = -1
      iloop = 0
      while(y<0 && iloop<nloop){
        x = iloop * step
        y = x*(1-x)**2-temp[i]
        iloop = iloop + 1
      }
      F1[i] = x
    }

    A_TAP = matrix(0, 1, N)
    for(i in 1:N){
      A_TAP[i] = A[i, i] * (1 - F1[i])
    }
    A_TAP = as.vector(A_TAP)
    A_TAP_inv = diag(A_TAP ^-1)

    W = A_TAP_inv %*% B
    W_thresh = threshold(W, threshold_type, ...)
    G = create_graph(W_thresh)

    structure(
      list(
        weights_matrix = W,
        thresholded_matrix = W_thresh,
        graph = G
      ),
      class = "ThoulessAndersonPalmer"
    )
}

cross_cov <- function(a, b){
  # cross_covariance
  # a,b -->  <(a - <a>)(b - <b>)>  (axis=0)

  da = sweep(a, 2, colMeans(a), '-')
  db = sweep(b, 2, colMeans(b), '-')
  r = (t(da) %*% db)/nrow(a)
  r
}

