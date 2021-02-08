#' # exact_mean_field.R
#' # ---------------------
#' #' Reconstruction of graphs using the exact mean field
#'
#' source(here::here( 'R', 'threshold.R'))
#'
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
#'   B = D %*% C_inv
#'
#'   if(exact) {
#'     f1 <- function(x, H) {
#'       (1 / sqrt(2*pi)) * exp(-(x**2)/2) * tanh(H+x*delta)
#'       }
#'     f2 <- function(x) {
#'       (1 / sqrt(2*pi)) * exp(-(x**2)/2) * (1 - ((tanh(H+x*delta))**2))
#'       }
#'
#'     W = matrix(NA, N, N)
#'     nloop = 100
#'     for (i0 in 1:N) {
#'       cost = matrix(0, 1, (nloop + 1))
#'       delta = 1
#'
#'       integrand <- function(H){
#'         y = stats::integrate(Vectorize(f1), -Inf, Inf, stop.on.error = F, H = H)$value
#'
#'         return(y - m[i0])
#'       }
#'
#'       for (iloop in 1:nloop) {
#'         #original python code is H = fsolve(integrand, 0.0)
#'         #H = stats::uniroot(integrand, c(-4,2))$root
#'         H = numDeriv::grad(integrand, x = 0)
#'         #H = nleqslv(0, integrand)
#'         #H = pracma::fsolve(integrand, 0)
#'
#'         a = stats::integrate(Vectorize(f2), -Inf, Inf, stop.on.error = F)$value
#'         print(a)
#'
#'
#'         if(a != 0){
#'           delta = (1/(a**2)) * sum(as.matrix(B[i0,]**2)*(1-m**2))
#'           W_temp = B[i0,] / a
#'           print("S")
#'           print(delta)
#'
#'           # sqrt(delta) will return nan if delta < 0
#'           # if(delta <  0){
#'           #   delta = delta
#'           # } else {
#'           #   delta = sqrt(delta)
#'           # }
#'
#'         } else {
#'           W_temp = B[i0,] *0
#'         }
#'
#'         H_temp = t(TS[, -L]) %*% W_temp
#'         cost[iloop] = mean((t(s1)[,i0]-tanh(H_temp))**2)
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
#'   #G = create_graph(W_thresh)
#'
#'   structure(
#'     list(
#'       weights_matrix = W,
#'       thresholded_matrix = W_thresh,
#'       graph = W_thresh
#'       #graph = G
#'     ),
#'     class = "MeanField"
#'   )
#' }
#'
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