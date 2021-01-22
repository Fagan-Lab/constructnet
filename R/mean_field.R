#' # exact_mean_field.R
#' # ---------------------
#' #' Reconstruction of graphs using the exact mean field
#' # author: Brennan Klein
#' # converted by: Zhaoyi Zhuang
#'
#' #' Title
#' #'
#' #' @param TS Array consisting of L observations from N sensors.
#' #' @param exact  If True, use the exact mean field approximation. If False, use the
#' #' naive mean field approximation.
#' #' @param stop_criterion  If True, prevent overly-long runtimes. Only applies for exact mean
#' #' field.
#' #' @param threshold_type Which thresholding function to use on the matrix of weights.
#' #' @param ... Arguments
#' #'
#' #' @return strcuture
#' #' @export
#' mean_field_fit <- function(TS, exact=T, stop_criterion=T, threshold_type='range', ...){
#'
#'   N = nrow(TS)
#'   L = ncol(TS)
#'   m = rowMeans(TS)
#'
#'   A = 1 - m ** 2
#'   A_inv = diag(1 / A)
#'   A = diag(A)
#'
#'   ds = sweep(t(TS), 2, m, '-')
#'   C = cov.wt(ds, method = "ML")$cov
#'   C_inv = solve(C)
#'
#'   s1 = TS[, -1]
#'   ds1 = sweep(t(s1), 2, rowMeans(s1), '-')
#'   D = cross_cov(ds1, ds[-L,])
#'
#'   # predict naive mean field W:
#'   B = D %*% C_inv
#'
#'   if(exact) {
#'     f1 <- function(x, H) (1 / sqrt(2*pi)) * exp(-(x**2)/2) * tanh(H+x*sqrt(delta))
#'     f2 <- function(x) (1 / sqrt(2*pi)) * exp(-(x**2)/2) * (1 - ((tanh(H+x*sqrt(delta)))**2))
#'
#'     W = matrix(NA, N, N)
#'     nloop = 100
#'     for (i0 in 1:N) {
#'       cost = matrix(0, 1, (nloop + 1))
#'       delta = 1
#'
#'       integrand <- function(H){
#'         y = pracma::quad(f1, -999, 999, H= H)
#'
#'         return(y - m[i0])
#'       }
#'
#'       for (iloop in 1:nloop) {
#'         H = integrand(1)
#'
#'         a = pracma::quad(f2, -999, 999)
#'
#'         if(a != 0){
#'           delta = (1/(a**2)) * rowSums(as.matrix(B[i0,]**2)*(1-m**2))
#'           W_temp = B[i0,] / a
#'         }
#'
#'         H_temp = t(TS[, -L]) %*% W_temp
#'         cost[iloop] = rowMeans((t(s1)[,i0]-tanh(H_temp))**2)
#'
#'         if(stop_criterion && iloop != 1 && (cost[iloop] >= cost[iloop - 1])){
#'           break
#'         }
#'       }
#'       W[i0, ] = W_temp
#'
#'     }
#'   } else {
#'     W  = A_inv %*% B
#'   }
#'   W_thresh = threshold(W, threshold_type, ...)
#'   G = create_graph(W_thresh)
#'
#'   structure(
#'     list(
#'       weights_matrix = W,
#'       thresholded_matrix = W_thresh,
#'       graph = G
#'     ),
#'     class = "MeanField"
#'   )
#'
#'
#' }
#'
#' cross_cov <- function(a, b) {
#'   # cross_covariance
#'   # a,b -->  <(a - <a>)(b - <b>)>  (axis=0)
#'
#'   da = sweep(a, 2, colMeans(a), '-')
#'   db = sweep(b, 2, colMeans(b), '-')
#'
#'   return((t(da) %*% db)/nrow(a))
#' }
#'
#'
#'
#'
#' TS = matrix(c(-1, -1, -2, -1, 4, 7, 8, -3, -2, 1, 1, 4, 0, 2, 2, 1, 3, 2, -3, -1, 9), nrow=3, ncol=7, byrow = TRUE)
#' TS
#' mean_field_fit(TS)
