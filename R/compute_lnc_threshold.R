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
#' @param lower Lower bound of constrained optimization
#' @param upper Upper bound of constrained optimization
#' @param num_iter Number of iterations of GP optimization
#' @param init_size Number of initial evaluation for estimating GP
#' @param diagnostics FALSE feature not implimented yet
#' @param cluster  A FORK cluster from parallel package
optimize_mse <- function(rho,
                         N,
                         M,
                         d,
                         k,
                         lower = -5,
                         upper = -1e-10,
                         num_iter = 25,
                         init_size = 20,
                         diagnostics = FALSE, #I would like to add some diagnostic features in the future
                         cluster = NULL) {

  objective_func <- function(alpha) {
    rmi::estimate_mse(k=k,alpha=alpha,d=d,rho=rho,N=N,M=M,cluster=cluster)
  }

  noise.var <- 0.1 #don't know what to do with this......

  # Generate initial exploration data
  alpha.doe   <- lower*as.data.frame(DiceDesign::lhsDesign(init_size, 1)$design)
  y.tilde     <- rep(0,init_size)

  cat("Evaluating initial design points...")

  for (i in 1:init_size) {
    y.tilde[i] <- objective_func(alpha.doe[[1]][i])
  }

  cat("Estimating kriging model...")

  # Estimate kriging model
  model <- DiceKriging::km(y~1, design=alpha.doe, response=data.frame(y=y.tilde),
              covtype="matern5_2", optim.method = "gen", noise.var=rep(noise.var,init_size),
              control=list(trace=FALSE))

  optim.param <- list()
  optim.param$quantile <- .9

  cat("Performing GP optimizations...")

  optim.result <- DiceOptim::noisy.optimizer(optim.crit="EQI", optim.param=optim.param, model=model,
                                  n.ite=num_iter, noise.var=noise.var, funnoise=objective_func,
                                  lower=lower, upper=upper,
                                  NoiseReEstimate=TRUE,
                                  CovReEstimate=TRUE,
                                  control=list(print.level = 0))

  best.x    <- optim.result$best.x

  return(best.x)
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
                         cluster = NULL) {

  inputs <- matrix(c(d,k,alpha,rho,N),ncol=5,nrow=M,byrow=TRUE)

  compute_mi <- function(input) {

    require(mnormt)
    simulate_mvn <- function(n,d,rho) { #could we utilize a seperate function for all of this?
      require(mnormt)
      Sigma       <- matrix(rho,d,d)
      diag(Sigma) <- 1
      return(rmnorm(n,mean(0,d),Sigma))
    }
    d = input[1]
    K = input[2]
    a = input[3]
    r = input[4]
    N = input[5]
    data   <- simulate_mvn(N,d,rho = r)
    return(rmi::knn_mi(data,splits = rep(1,d), options = list(method="LNC",k=K,alpha=c(a,rep(0,d)))))
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

  return(mean((mi_mse_est - analytic_mi(d,rho))^2))
}










