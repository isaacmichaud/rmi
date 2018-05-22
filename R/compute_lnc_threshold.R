#' Title
#'
#' @param d
#' @param k
#' @param cached
#'
#' @return
#' @export
#'
#' @useDynLib rmi
#' @importFrom Rcpp sourceCpp
#'
#' @examples
compute_lnc_threshold <- function(d,k,cached=TRUE) {
 if ((d <= 100) & (k <= 100)) {
    return(rmi::alpha_thresholds[d,k])
 } else {
  #brute force that stuff.....
 }
}
