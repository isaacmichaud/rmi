#' Random Multivariate Normal (MVN)
#'
#' @param n sample size
#' @param mu mean vector
#' @param cov_mat covariance matrix
#'
#' @return Matrix (n x d) of random MVN draws
#'
#' @examples
#' mnormt::rmnorm(200,0, matrix(c(1,0.9,0.9,1),2,2))
rmvn <- function(n,mu = rep(0,d),cov_mat) {
  S <- chol(cov_mat)
  d <- ncol(S)
  Y <- (matrix(rnorm(n*d),n,d) %*% S) + mu
  return(Y)
}

#' Calibration Random Multivariate Normal
#'
#' @param n integer, sample size
#' @param d integer, dimension
#' @param rho correlation of compound-symmetric covariance
#'
#' @return Matrix (n x d) of random MVN draws
#'
#' @examples
#' simulat_mvn(100,10,0.6)
simulate_mvn <- function(n,d,rho) {
  Sigma       <- matrix(rho,d,d)
  diag(Sigma) <- 1
  return(rmi::rmvn(n,mean(0,d),Sigma))
}
