#' LNN Entropy Estimator
#'
#' Local Nearest Neighbor entropy estimator using Gaussian kernel and kNN selected bandwidth. Entropy is estimated by taking
#'
#' @param data matrix, each row is an observation
#' @param k  order of the local kNN bandwidth selection
#' @param tr number of neighbors to include in local density estimate (truncation)
#' @param bw bandwidth (optional) manually fix bandwidth instead of using local kNN bandwidth selection
#'
#' @section References:
#'
#' Loader, C. (1999). Local regression and likelihood. Springer Science & Business Media.
#'
#' Gao, W., Oh, S., & Viswanath, P. (2017). Density functional estimators with k-nearest neighbor bandwidths. IEEE International Symposium on Information Theory - Proceedings, 1, 1351–1355. https://doi.org/10.1109/ISIT.2017.8006749
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(1000)
#' print(lnn_entropy(x))
#' #analytic entropy
#' print(0.5*log(2*pi*exp(1)))
#'
#' @export
lnn_entropy <- function(data, k = 5, tr = 30, bw=NULL) {

    #bandwidth selection is k, number of points used in the calculation is tr
    if (k > tr) {
      stop("Unaccaptable Inuputs, k is larger than tr")
    }

    if (!is.matrix(data)) {
      data <- matrix(data,ncol=1)
    }

    nn <- nearest_neighbors(data,tr)
    N  <- dim(data)[1]
    d  <- dim(data)[2]

    local_estimate <- rep(0,N)

    if(is.null(bw)) {
      bw = nn$nn_dist[,k+1]
    } else {
      bw = rep(bw,N)
    }

    for (i in 1:N) {
      S0 = 0
      S1 = rep(0,d)
      S2 = matrix(0,d,d)
      for (j in 1:(tr+1)) {
        dis = t(data[nn$nn_inds[i,j],] - data[i,])
        w   = as.numeric(exp(-(dis%*%t(dis))/(2*bw[i]^2)))
        S0 = S0 + w
        S1 = S1 + w*(dis/bw[i])
        S2 = S2 + w*((t(dis)%*%dis)/(bw[i]^2))
      }
      Sigma = S2/S0 - (t(S1)%*%S1)/(S0^2)

      if (is.nan(Sigma)) {
        stop("Local covariance matrix was singular")
      }

      if (d == 1) {
        det_Sigma = as.numeric(Sigma)
      } else {
        det_Sigma = det(Sigma)
      }

      if (det_Sigma < (1e-4)^d) {
        local_estimate[i] = 0
      } else {
        if (d == 1) {
          offset = (S1/S0)%*%(1/as.numeric(Sigma)*(S1/S0))
        } else {
          offset = t((S1/S0))%*%t(solve(Sigma)%*%t((S1/S0)))
        }
        local_estimate[i] = -log(S0) + log(N-1) + 0.5*d*log(2*pi) + d*log(bw[i]) + 0.5*log(det_Sigma) + 0.5*offset[1]
      }
    }

    if (sum(local_estimate > 0) == 0) {
      return(0)
    } else {
      return(mean(local_estimate[local_estimate != 0]))
    }
    warning("Bad returned value")
    return(NA)
}

#' LNN MI estimator
#'
#' Local nearest neighbor mutual information estimator using Gaussian kernel and kNN bandwidth selection. Plug-in estimate of mutual information using LNN entropy estimator.
#'
#' @param data
#' @param splits
#' @param k  order of the local kNN bandwidth selection
#' @param tr number of neighbors to include in local density estimate (truncation)
#'
#' @export
lnn_mi <- function(data, splits, k = 5, tr = 30) {
  xy_E <- lnn_entropy(cbind(x,y),k=k,tr=tr)
   x_E <- lnn_entropy(x,k=k,tr=tr)
   y_E <- lnn_entropy(y,k=k,tr=tr)
   return(x_E + y_E - xy_E)
}
