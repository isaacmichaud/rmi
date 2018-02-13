library(rmi)
context("LNC Estimates")

test_that("LNC does something",{
  set.seed(123)
  n = 5000
  T = matrix(c(1,0.7,0.7,1),ncol=2,nrow=2)
  x = t(T%*%matrix(rnorm(2*n),nrow=2))
  mi = knn_mi(x,c(1,1),options = list(method = "LNC", k = 10, alpha = c(0.6868716,0,0)))
  TT = t(T)%*%T
  r  = TT[1,2]/TT[2,2]
  expect_equal(mi,-0.5*log(1-r^2))
}
)


simulate_mvn <- function(n,d,rho) {
  require(mnormt)
  Sigma       <- matrix(rho,d,d)
  diag(Sigma) <- 1
  return(rmnorm(n,mean(0,d),Sigma))
}

analytic_mi <- function(d,rho,splits=c(2,2)) {
  Sigma       <- matrix(rho,d,d)
  diag(Sigma) <- rep(1,d)
  S22_inv     <- solve(Sigma[-(1:splits[1]),-(1:splits[1])])
  S12         <- Sigma[(1:splits[1]),-(1:splits[1])]
  S21         <- Sigma[-(1:splits[1]),(1:splits[1])]
  cond_Sigma  <- Sigma[1:splits[1],1:splits[1]] - S12%*%S22_inv%*%S21
  return(0.5*log(det(Sigma[(1:splits[1]),(1:splits[1])])/det(cond_Sigma)))
}

expect_equal(analytic_mi(2,0.9,c(1,1)),-0.5*log(1-0.9^2))

set.seed(123)
n = 5000
x = simulate_mvn(n,10,0.99)
mi = knn_mi(x,c(7,3),options = list(method = "LNC", k = 15, alpha = c(0.01272956,0,0)))
expect_equal(mi,analytic_mi(10,0.99,splits=c(3,7)))
