compute_lnc_threshold <- function(d,k,cached=TRUE) {
 if ((d <= 100) & (k <= 100)) {
    return(rmi::alpha_thresholds[d,k])
 } else {
  #brute force that stuff.....
 }
}
