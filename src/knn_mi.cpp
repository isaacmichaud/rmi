#include <RcppArmadillo.h>
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
#include "parse_split_vector.h"
#include "get_nearest_neighbors.h"
#include "lnc_compute.h"
using namespace std;

//' mutual information function
//'
//' @param x numeric vector
//' @param low lower bound
//' @param high upper bound
//' @param res bounded numeric vector
//' @return bounded numeric vector
//' @export
//'
// [[Rcpp::export]]
double knn_mi(arma::mat data,
              Rcpp::NumericVector splits,
              const Rcpp::List & options) {

  std::string method = Rcpp::as<std::string>(options["method"]);
  int              k = Rcpp::as<int>(options["k"]);
  int            lnc = 0;

  int K     = k+1;
  int N     = data.n_rows;
  int vars  = splits.length();
  arma::colvec alpha(vars+1);

  if (method == "LNC") {
    lnc    = 1;
    method = "KSG2";
    alpha  = Rcpp::as<arma::colvec>(options["alpha"]);
  }

  arma::imat nn_inds(N,K);
  arma::mat nn_dist(N,K);

  arma::icolvec d(vars+1);
  arma::icolvec d_start(vars+1);
  arma::icolvec d_end(vars+1);

  parse_split_vector(splits,d,d_start,d_end);
  get_nearest_neighbors(data, nn_dist, nn_inds,k);
  double digamma_x = 0;
  double mi = (vars-1)*boost::math::digamma(N) + boost::math::digamma(k);
  int N_x;
  double epsilon;
  double dist;
  double proposed_correction;
  double lnc_correction = 0;

  if (method == "KSG1") { //ksg1
    for (int i = 0; i < N; i++) { // for point i
      for (int j = 1; j <= vars; j++) { // count marginal sum for jth block
        N_x     = 0;
        epsilon = nn_dist(i,k);

        for (int m = 0; m < N; m++) { // iterate over all points
          dist = 0;

          for (int n = d_start(j); n <= d_end(j); n++) {
            if (dist < abs(data(i,n)-data(m,n))) {
              dist = abs(data(i,n)-data(m,n));
            }
          }
          if (dist < epsilon) N_x++;
        }
        digamma_x = digamma_x + boost::math::digamma(N_x+1);
      }
    }
   // printf("mi: %f, digamma:%f",mi,digamma_x);
    mi = mi - digamma_x/(double)N;
  }

  if (method == "KSG2") { //ksg2
    for (int i = 0; i < N; i++) { // for point i
      for (int j = 1; j <= vars; j++) { // count marginal sum for jth block
        N_x = 0;
        epsilon = 0;

        for (int m = 1; m < K; m++) {
          dist = 0;
          for (int n = d_start(j); n <= d_end(j); n++) {
            if (dist < abs(data(i,n)-data(nn_inds(i,m),n))) {
              dist = abs(data(i,n)-data(nn_inds(i,m),n));
            }
          }
          if (dist > epsilon) {
            epsilon = dist;
          }
        }

        for (int m = 0; m < N; m++) { // iterate over all points
          // printf("%d\n",j);
          dist = 0;

          for (int n = d_start(j); n <= d_end(j); n++) {
            if (dist < abs(data(i,n)-data(m,n))) {
              dist = abs(data(i,n)-data(m,n));
            }
          }
          //  printf("dist:%f\n",dist);
          if (dist <= epsilon) N_x++;
        }
        digamma_x = digamma_x + boost::math::digamma(N_x);
      }
    }
    //printf("mi: %f, digamma:%f",mi,digamma_x);
    mi = mi - digamma_x/(double)N - (vars - 1)/(double)k ;
  }

  // if (lnc == 1) { //LNC corrections
  //   for (int j = 0; j < N; j++) {
  //     proposed_correction = lnc_compute(XY, XY_inds, i,d,K);
  //     if (proposed_correction < log(alpha(0))) {
  //       // printf("%d : %f, %f\n",j,proposed_correction,log(alpha(0)));
  //       lnc_correction = lnc_correction - proposed_correction;
  //     }
  //     for (int i = 1; i <= vars; i++) {
  //       if (splits(i-1) == 1) continue;
  //       proposed_correction = lnc_compute(X, X_inds, i,d_x,K);
  //       if (proposed_correction < alpha(i)) {
  //         // printf("%d : %f, %f\n",j,proposed_correction,log(alpha(i)));
  //         lnc_correction = lnc_correction - proposed_correction;
  //       }
  //     }
  //   }
  // }

  if (lnc == 1) { //LNC corrections
    for (int i = 0; i < N; i++) {
      for (int j = 0; j <= vars; j++) {
        //skip any coordinates with 1 variable
        if (d_end(j) - d_start(j) == 0) continue;
        proposed_correction = lnc_compute(data, nn_inds, i, d_start(j), d_end(j));
        if (proposed_correction < log(alpha(j))) {
          // printf("%d : %f, %f\n",j,proposed_correction,log(alpha(0)));
          if (j == 0) {
            lnc_correction = lnc_correction - proposed_correction;
          } else {
            lnc_correction = lnc_correction + proposed_correction;
          }
        }
      }
    }
  }

  return mi + lnc_correction/(double)N;
}

