library(rmi)
context("kNN estimates")

test_that("knn does something",{
  set.seed(123)
  x = rnorm(1000)
  y = rnorm(1000) + x
  mi = knn_mi(cbind(x,y),c(1,1),options = list(method = "KSG2", k = 5))
  expect_equal(mi,0.349905528753732)
}
)
