#' Title
#'
#' @param d
#' @param k
#' @param cached
#'
#' @return x
#' @export
#' @examples x+y = 2
compute_lnc_threshold <- function(d,k,cached=TRUE) {
 if ((d <= 100) & (k <= 100)) {
    return(rmi::alpha_thresholds[d,k])
 } else {
  #brute force that stuff.....
 }
}
