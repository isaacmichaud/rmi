#include <RcppArmadillo.h>
#include "flann/flann.hpp"

using namespace std;

void get_nearest_neighbors(arma::mat, arma::mat&, arma::imat&, int);

Rcpp::List nearest_neighbors(arma::mat, int);
