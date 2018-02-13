#include <RcppArmadillo.h>
#include "flann/flann.hpp"
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
using namespace std;

//' test knn algorithm
//'
//' @param x numeric vector
//' @param low lower bound
//' @param high upper bound
//' @param res bounded numeric vector
//' @return bounded numeric vector
//' @export
// [[Rcpp::export]]
Rcpp::List test_knn(Rcpp::NumericMatrix data, int N, int k) {
  int K             = k+1;
  arma::mat nn_dist = arma::zeros(N,K);
  arma::mat nn_inds = arma::zeros(N,K);
  const int d       = data.ncol();
  arma::mat query(d, N);
  arma::mat temp_q(data.begin(), N, d, false);
  query = arma::trans(temp_q); //transpose
  flann::Matrix<double> q_flann(query.memptr(), N, d);
  arma::mat ref(d, N);
  {
    arma::mat temp_r(data.begin(), N, d, false);
    ref = arma::trans(temp_r);
  }
 flann::Matrix<double> ref_flann(ref.memptr(), N, d);
 // // Setting the flann index params
 flann::IndexParams params;
 // if (build == "kdtree") {
 //params = flann::KDTreeSingleIndexParams(1);
 // } else if (build == "kmeans") {
 //   params = flann::KMeansIndexParams(2, 10, flann::FLANN_CENTERS_RANDOM, 0.2);
 // } else if (build == "linear") {
 params = flann::LinearIndexParams();
 //params = flann::KMeansIndexParams();
 // }
 // // Finding the nearest neighbours
 //flann::Index<flann::L2<double> > index(ref_flann, params);
 flann::Index<flann::MaxDistance<double> > index(ref_flann, params);
 index.buildIndex();
 flann::Matrix<int> indices_flann(new int[N * k], N, k);
 flann::Matrix<double> dists_flann(new double[N * k], N, k);
 flann::SearchParams search_params;
 search_params.cores = 1; //change later
 search_params.checks = 100;
 index.knnSearch(q_flann, indices_flann, dists_flann, k, search_params);
 arma::imat indices(indices_flann.ptr(), k, N, true);
 arma::mat dists(dists_flann.ptr(), k, N, true);
 delete[] indices_flann.ptr();
 delete[] dists_flann.ptr();
 return Rcpp::List::create(Rcpp::Named("indices")   = indices.t() + 1,
                           Rcpp::Named("distances") = dists.t());
 // return Rcpp::List::create(Rcpp::Named("nn_dist") = 0,Rcpp::Named("nn_inds") = 0);
  //return Rcpp::List::create(Rcpp::Named("nn_dist") = nn_dist,Rcpp::Named("nn_inds") = nn_inds);
}



