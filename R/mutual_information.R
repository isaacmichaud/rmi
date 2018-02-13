#' k Nearest Neighbor Mutual Information Estimation
#'
#' Here is an explaination of what the stuff does
#'
#' @param data numeric vector
#' @param splits lower bound
#' @param k upper bound
#' @param method bounded numeric vector
#' @return bounded numeric vector
#' @export
#'
knn_mi_temp <- function(data, splits, k, method = "KSG2") {
  #--- checks ---#
  N <- length(data[1,])
  if (k > N) stop("k is too large for the supplied number of data points")
  if (k > (N/2)) warning("k is larger than half of supplied number of data")

  d <- length(data[,1]) #joint dimension
  if (sum(splits) != joint_dimension) stop("Data dimensionality is not the same as implied by splits vector")

  #--- initializations ---#
  if (method == "LNC") {
    alpha    <- rep(0,length(splits)+1)
    alpha[1] <- rmi::get_lnc_thresholds(d,k)
    for (i in seq_along(splits)) {
      alpha[i+1] <- rmi::get_lnc_thresholds(d,k)
    }
    options = list(method = method, k = k, alpha = alpha)
  } else {
    options = list(method = method, k = k)
  }
  #--- call c++ code---#
  .Call('_rmi_knn_mi', PACKAGE = 'rmi', data, splits, options)
}

