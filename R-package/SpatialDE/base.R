#' Squared Exponential Kernel
#'
#' Gaussian covariance function for the Gaussian Processs model, for general
#' spatial covariance.
#'
#' @param X Gene expression matrix with genes in rows and cells in columns.
#' @param l Controls how wide the kernel is; larger \code{l} means wider.
#' @param threads Number of CPU threads for calculating the disstance matrix. Default
#' is the maximum number of CPU threads available on the system.
#' @family kernels
#' @return A covariance matrix for cells in \code{X}.
#' @seealso \code{\link{linear_kernel}}, \code{\link{cosine_kernel}}

SE_kernel <- function(X, l, threads = NULL) {
  R2 <- as.matrix(parallelDist::parDist(X, threads = threads)) ^ 2
  exp(-R2 / (2 * l ^ 2))
}

#' Linear kernel
#'
#' Linear covariance function for the Gaussian Processs model.
#'
#' @param X Gene expression matrix with genes in rows and cells in columns.
#' @param threads Number of CPU threads for calculating the disstance matrix. Default
#' is the maximum number of CPU threads available on the system.
#' @family kernels
#' @return A covariance matrix for cells in \code{X}.
#' @seealso \code{\link{SE_kernel}}, \code{\link{cosine_kernel}}
#'
linear_kernel <- function(X) {
  X %*% t(X)
}

#' Cosine kernel
#'
#' Cosine covariance function for the Gaussian Processs model, for periodic spatial
#' patterns.
#'
#' @param X Gene expression matrix with genes in rows and cells in columns.
#' @param p Period of the cosine.
#' @param threads Number of CPU threads for calculating the disstance matrix. Default
#' is the maximum number of CPU threads available on the system.
#' @family kernels
#' @return A covariance matrix for cells in \code{X}.
#' @seealso \code{\link{SE_kernel}}, \code{\link{linear_kernel}}
#'
cosine_kernel <- function(X, p) {
  dist_mat <- as.matrix(parallelDist::parDist(X, threads = threads))
  cos(pi/p * dist_mat)
}

#' Gower Scaling Factor
#'
#' Gower normalization factor for covariance matrix \code{K}, used when estimating the
#' proportion of variance explained by spatial variance.
#'
#' @param K Covariance matrix.
#' @return The Gower scaling factor
#'
gower_scaling_factor <- function(K) {
  n <- nrow(K)
  P <- diag(nrow = n) - 1/n
  sum(diag(P %*% K %*% P)) / (n - 1)
}

factor <- function(K) {
  factors <- eigen(K)
  factors$values[factors$values < 0.] <- 0.
  list('U' = factors$vectors, 'S' = factors$values)
}

get_UT1 <- function(U) {
  colSums(U)
}

get_UTy <- function(U, y) {
  t(y) %*% U
}

mu_hat <- function(delta, UTy, UT1, S) {
  UT1_scaled <- UTy / (S + delta)
  sum1 <- UT1_scaled %*% t(UTy)
  sum2 <- UT1_scaled %*% UT1
  sum1 / sum2
}

LL <- function(delta, UTy, UT1, S, n) {
  mu_h <- mu_hat(delta, UTy, UT1, S)

  sum_1 <- sum((UTy - UT1 * mu_h) / (S + delta))
  sum_2 <- sum(log(S + delta))

  -0.5 * (n * log(2 * pi) + n * log(sum_1 / n) + sum_2 + n)
}

make_objective <- function(UTy, UT1, S, n){
  ll_obj <- function(log_delta){
    -LL(exp(log_delta), UTy, UT1, S, n)
  }
  ll_obj
}

lengthscale_fit <- function(UTy, UT1, S, n) {
  ll_obj <- make_objective(UTy, UT1, S, n)
  o <- optimize(ll_obj, c(-10, 20))
  list('max_delta' = exp(o$minimum), 'max_ll' = -o$objective)
}

lengthscale_fit_helper <- function(y, U, UT1, S) {
  UTy <- get_UTy(U, y)
  n <- length(y)
  lengthscale_fit(UTy, UT1, S, n)
}
