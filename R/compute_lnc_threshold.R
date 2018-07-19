#' Compute LNC Threshold
#'
#' Computes alpha tuning parameters using a fitted surrogate model. Uses the minimize MSE method of calulcating alpha
#'
#' @param d number of dimenions
#' @param k number of neighbors
#' @param n sample size
#'
#' @return log-scale value of alpha value will be on the interval (-inf,0)
#' @export
#'
#' @examples
#' compute_lnc_threshold(2,5,1000)
#'
compute_lnc_threshold <- function(d, k, n) {
  warning("surrogate threshold computer is not implimented yet")
  return(NA)
}

#' Optimize MSE of alpha MI estimators
#' Uses Gaussian process (GP) optimization to stochastically minimize the mean squared error of the LNC estimator
#'
#' @param rho Correlation of reference distribution
#' @param N   Sample size
#' @param M   Monte Carlo sample size
#' @param d   Dimension
#' @param k   Neighboorhood order
#' @param num_iter Number of iterations of GP optimization
#' @param init_size Number of initial evaluation for estimating GP
#' @param cluster  A FORK cluster from parallel package
#' @export
optimize_mse <- function(rho,
                         N,
                         M,
                         d,
                         k,
                         lower = -10,
                         upper = -1e-10,
                         num_iter = 10,
                         init_size = 20,
                         diagnostics = FALSE, #I would like to add some diagnostic features in the future
                         cluster = NULL,
                         verbose = TRUE) {

  objective_func <- function(alpha) {
    estimate_mse(k=k,alpha=alpha,d=d,rho=rho,N=N,M=M,save_result = TRUE)
  }

  #--- Find Transition Point ---#
  rect      <- rbind(c(lower,upper))
  out       <- NULL
  progress  <- NULL
  threshold <- NULL
  X         <- NULL
  Z         <- NULL
  while (is.null(threshold)) {
    Xcand  <- tgp::lhs(init_size, rect)
    Xnew  <- tgp::dopt.gp(init_size, X = X, Xcand)$XX
    X     <- rbind(X, Xnew)
    new_Z <- NULL
    for (i in seq_along(Xnew)) {
      new_Z[i] <- objective_func(Xnew[i,])
    }
    Z   <- c(Z, new_Z)
    out <- tgp::optim.step.tgp(objective_func, X=X, Z=Z, rect=rect, prev=out)
    threshold <- tryCatch({out$obj$trees[[2]]$val[1]},error=function(e){NULL})
    if (min(Z) > 0.5) {
      rect[1] <- 2*rect[1]
    }
  }

  #--- Improve Design ---#
  rect <- threshold + c(-1,1)
  Xcand  <- tgp::lhs(init_size, rect)
  Xnew  <- tgp::dopt.gp(init_size, X = X, Xcand)$XX
  X     <- rbind(X, Xnew)
  new_Z <- NULL
  for (i in seq_along(Xnew)) {
    new_Z[i] <- objective_func(Xnew[i,])
  }
  Z   <- c(Z, new_Z)

  #--- Adaptive Optimization ---#
  for(j in 1:num_iter) {

    out <- tgp::optim.step.tgp(objective_func, X=X, Z=Z, rect=rect, prev=out)

    ## add in the inputs, and newly sampled outputs
    X <- rbind(X, out$X)
    new_Z <- NULL
    for (i in seq_along(out$X)) {
      new_Z[i] <- objective_func(out$X[i,])
    }
    Z <- c(Z, new_Z)

    ## keep track of progress and best optimum
    progress <- rbind(progress, out$progress)
    if (verbose) {
      print(paste(sprintf("Iteration %d of %d :",j,num_iter),print(out$progress$x1),sep = ""))
    }
  }
  if (verbose) plot(out$obj)
  return(out$progress$x1)
}


#' Estimate MSE of KNN Estimator
#'
#' This function computes the MSE for a given set of tuning parameter alpha, dimension, size of neighborhood
#'
#' @param k
#' @param alpha alpha if < 0 then log scale otherwise [0,1]
#' @param d dimension
#' @param rho MVN correlation (compound systemtric)
#' @param N sample size
#' @param M replications (outer loop)
#' @param cluster a parallel cluster object
#'
#' @return
#' @export
#'
#' @examples
estimate_mse <- function(k       = 5,
                         alpha   = 0,
                         d       = 2,
                         rho     = 0.0,
                         N       = 1000,
                         M       = 100,
                         cluster = NULL,
                         save_result = FALSE) {

  inputs <- matrix(c(d,k,alpha,rho,N),ncol=5,nrow=M,byrow=TRUE)

  compute_mi <- function(input) {

    simulate_mvn <- function(n,d,rho) {
      Sigma       <- matrix(rho,d,d)
      diag(Sigma) <- 1
      return(rmvn(n,mean(0,d),Sigma))
    }
    d = input[1]
    K = input[2]
    a = input[3]
    r = input[4]
    N = input[5]
    data   <- simulate_mvn(N,d,rho = r)
    return(knn_mi(data,splits = rep(1,d), options = list(method="LNC",k=K,alpha=c(a,rep(0,d)))))
  }

  if (is.null(cluster)) {
    mi_mse_est <- rep(0,M)
    for (i in 1:M) {
      mi_mse_est[i] <- compute_mi(inputs[i,])
    }
  } else {
    mi_mse_est <- parApply(cluster,inputs,1,compute_mi)
  }

  analytic_mi <- function(d,rho) { #this would be a good function to break off too
    Sigma       <- matrix(rho,d,d)
    diag(Sigma) <- 1
    return(-0.5*log(det(Sigma)))
  }
  my_mse <- mean((mi_mse_est - analytic_mi(d,rho))^2)
  if (save_result) {
    save_val <- c(k, alpha, d, rho, N, M, my_mse)
    write(save_val,file = "lnc_alpha_archive.txt",append = TRUE)
  }
  return(my_mse)
}











