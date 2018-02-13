#' mutual information test
#'
#' @param x numeric vector
#' @param low lower bound
#' @param high upper bound
#' @param res bounded numeric vector
#' @return bounded numeric vector
#' @export
run_mi_test = function() {
  set.seed(1234)
  res = rep(0,1000)
  for (i in 1:1000) {
    X = rnorm(1000)
    Y = rnorm(1000) + X
    res[i] = rmi::mutual_information(as.matrix(X),as.matrix(Y),cbind(X,Y),5,2)
  }
  print(summary(res))
  hist(res)
}

#' test centering
#'
#' @param x numeric vector
#' @param low lower bound
#' @param high upper bound
#' @param res bounded numeric vector
#' @return bounded numeric vector
#' @export
run_center_test = function() {
  set.seed(1234)
  X = rnorm(1000)
  Y = rnorm(1000) + X
  temp = test_get_nearest_neighbors(as.matrix(cbind(X,Y)),5)
  res = compute_LNC(as.matrix(cbind(X,Y)),temp[[2]],0,2,6)
  plot(X,Y)
  points(X[temp[[2]][1,]+1],Y[temp[[2]][1,]+1],col='red',pch=20)
  points(res[,1],res[,2],col='green')
}
