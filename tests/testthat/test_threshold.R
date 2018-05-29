library(rmi)
context("Threshold Estimation")

test_that("Parallel execution works",{
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

test_that("log_alpha values give same estimates",{
  n      <- 5000
  T      <- matrix(c(1,0.7,0.7,1),ncol=2,nrow=2)
  x      <- t(T%*%matrix(rnorm(2*n),nrow=2))
  mi     <- knn_mi(x,c(1,1),options = list(method = "LNC", k = 10, alpha = c(0.6868716,0,0)))
  mi_log <- knn_mi(x,c(1,1),options = list(method = "LNC", k = 10, alpha = c(log(0.6868716),0,0)))
  expect_equal(mi,mi_log)
}
)
