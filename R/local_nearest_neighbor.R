#' kNN Mutual Information Estimators

#' @param  data matrix. Each row is an observation.
#' @param  splits vector. Describes which sets of columns in \code{data} to compute the mutual information between. For example, to compute mutual information between two variables use \code{splits = c(1,1)}. To compute \emph{redundancy} among multiple random variables use \code{splits = rep(1,ncol(data))}. To compute the mutual information between two random vector list the dimensions of each vector.
#' @param  options list. Specifies estimator and necessary parameters. See Details.
#' @section Details: Current types of methods that are available are methods LNC, KSG1 and KSG2
#' \code{list(method = "KSG2", k = 6)}
#' @return estimated mutual information
#' @section Author:
#' Isaac Michaud, North Carolina State University, \email{ijmichau@ncsu.edu}
#' @section References:
#' Gao, Shuyang, Greg Ver Steeg, and Aram Galstyan. 2015. "Efficient estimation of mutual information for strongly dependent variables." Artificial Intelligence and Statistics: 277-286.
#'
#' Kraskov, Alexander, Harald Stogbauer, and Peter Grassberger. 2004. "Estimating mutual information." Physical review E 69(6): 066138.
#'


#' LNN Entropy Estimator
#'
#' Computes mutual information based on the distribution of nearest neighborhood distances as described by Kraskov, et. al (2004).
#'
#' @param data
#' @param k
#' @param tr
#' @param bw
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(1000)
#' y <- x + rnorm(1000)
#' knn_mi(cbind(x,y),c(1,1),options = list(method = "KSG2", k = 6))
#'
#' set.seed(123)
#' x <- rnorm(1000)
#' y <- 100*x + rnorm(1000)
#' knn_mi(cbind(x,y),c(1,1),options = list(method = "LNC", alpha = 0.65, k = 10))
#' #approximate analytic value of mutual information
#' -0.5*log(1-cor(x,y)^2)
#'
#' z <- rnorm(1000)
#' #redundancy I(x;y;z) is approximately the same as I(x;y)
#' knn_mi(cbind(x,y,z),c(1,1,1),options = list(method = "LNC", alpha = c(0.5,0,0,0), k = 10))
#' #mutual information I((x,y);z) is approximately 0
#' knn_mi(cbind(x,y,z),c(2,1),options = list(method = "LNC", alpha = c(0.5,0.65,0), k = 10))
#'
#'
#' @export
lnn_entropy <- function(data,k=5,tr=30,bw=NULL) {
    #bandwidth selection is k, number of points used in the calculation is tr
    if (k > tr) {
      stop("Unaccaptable Inuputs, k is larger than tr")
    }

    nn <- rmi::nearest_neighbors(data,tr)
    N  <- dim(data)[1]
    d  <- dim(data)[2]

    local_estimate <- rep(0,N)

    if(is.null(bw)) {
      bw = nn$nn_dist[,k+1]
    } else {
      bw = rep(bw,N)
    }

    for (i in 1:N) {
      #lists = tree.query(x[i],tr+1,p=2)
      #knn_dis = lists[0][k]
      S0 = 0
      S1 = rep(0,d)
      S2 = matrix(0,d,d)
      #browser(expr = (i == 52))
      for (j in 1:(tr+1)) {
        #browser(j == 10)
        dis = t(data[nn$nn_inds[i,j],] - data[i,])
        w   = as.numeric(exp(-(dis%*%t(dis))/(2*bw[i]^2)))
        S0 = S0 + w
        S1 = S1 + w*(dis/bw[i])
        S2 = S2 + w*((t(dis)%*%dis)/(bw[i]^2))
      }
      #browser()
      Sigma = S2/S0 - (t(S1)%*%S1)/(S0^2)
      #print(c(Sigma,w))
      browser(expr = is.nan(Sigma))

      if (d == 1) {
        det_Sigma = as.numeric(Sigma)
      } else {
        det_Sigma = det(Sigma)
      }
      #print(c(det_Sigma,(1e-4)^d))
      if (det_Sigma < (1e-4)^d) {
        local_estimate[i] = 0
      } else {
        if (d == 1) {
          offset = (S1/S0)%*%(1/as.numeric(Sigma)*(S1/S0))
        } else {
          offset = t((S1/S0))%*%t(solve(Sigma)%*%t((S1/S0)))
        }
        local_estimate[i] = -log(S0) + log(N-1) + 0.5*d*log(2*pi) + d*log(bw[i]) + 0.5*log(det_Sigma) + 0.5*offset[1]
        #local_estimate[i] = -log(S0) + log(N-1) + 0.5*d*log(2*pi) + d*log(bw[i]) + 0.5*log(det_Sigma)
      }
    }
    #browser()
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
#' @param x
#' @param y
#' @param k
#' @param tr
#'
#' @export
#'
lnn_mi <- function(x,y,k = 5, tr = 30) {
  xy_E <- lnn_entropy(cbind(x,y),k=k,tr=tr)
   x_E <- lnn_entropy(x,k=k,tr=tr)
   y_E <- lnn_entropy(y,k=k,tr=tr)
   return(x_E + y_E - xy_E)
}

#I don't know if the bias term is correctly calibrated in this code. The MI estimate is close, but not really that good
#x = matrix(rnorm(1000),ncol=1)
#y = x + rnorm(1000)
#print(lnn_mi(x,y,tr=100))


