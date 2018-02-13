#include <RcppArmadillo.h>
#include "flann/flann.hpp"
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
using namespace std;


void get_nearest_neighbors(Rcpp::NumericMatrix X,
                           arma::mat& X_dist,
                           arma::imat& X_inds,
                           int k) {
  int d = X.ncol();
  int N = X.nrow();
  int K = k + 1;

  arma::mat query(d, N);
  arma::mat temp_q(X.begin(), N, d, false);
  query = arma::trans(temp_q); //transpose
  flann::Matrix<double> q_flann(query.memptr(), N, d);
  flann::Matrix<double> ref_flann(query.memptr(), N, d);

  //Setting the flann index params
  flann::IndexParams params;
  params = flann::LinearIndexParams();

  flann::Index<flann::MaxDistance<double> > index(ref_flann, params);
  index.buildIndex();
  flann::Matrix<int> indices_flann(new int[N * K], N, K);
  flann::Matrix<double> dists_flann(new double[N * K], N, K);
  flann::SearchParams search_params;

  //search_params.cores = 10;
  //search_params.checks = 100;

  index.knnSearch(q_flann, indices_flann, dists_flann, K, search_params);

  arma::imat indices(indices_flann.ptr(), K, N, true);
  arma::mat dists(dists_flann.ptr(), K, N, true);

  delete[] indices_flann.ptr();
  delete[] dists_flann.ptr();

  X_inds = indices.t();
  X_dist = dists.t();

  return;
}

//' test get nearest function
//'
//' @param x numeric vector
//' @param low lower bound
//' @param high upper bound
//' @param res bounded numeric vector
//' @return bounded numeric vector
//' @export
// [[Rcpp::export]]
double compute_LNC(Rcpp::NumericMatrix X, arma::imat X_inds, int index, int d, int K) {
  // initial declarations
  arma::colvec mean = arma::zeros(d);
  arma::mat data = arma::zeros(K,d);
  double correction = 0;

  // center the data
  for (int i = 0; i < K ; i++) {
    for (int j = 0; j < d ; j++) {
      mean(j) = mean(j) + X(X_inds(index,i),j)/K;
    }
  }

  for (int i = 0; i < K ; i++) {
    for (int j = 0; j < d ; j++) {
      data(i,j) = X(X_inds(index,i),j) - mean(j);
    }
  }

  // compute coordinate aligned volume
  double rect_volume = 0;
  arma::colvec rect_dist = arma::zeros(d);
  for (int i = 0; i < K ; i++) {
    for (int j = 0; j < d ; j++) {
      if (rect_dist(j) < abs(data(i,j))) {
        rect_dist(j) = abs(data(i,j));
      }
    }
  }

  for (int j = 0; j < d ; j++) {
    rect_volume = rect_volume + log(rect_dist(j));
  }

  // compute PCA aligned volume
  arma::mat U(K,d);
  arma::mat V(d,d);
  arma::colvec s(d);
  arma::svd(U, s, V, data);

  arma::mat pca_data(K,d);
  pca_data = data*V;

  double pca_volume = 0;
  arma::colvec pca_dist = arma::zeros(d);
  for (int i = 0; i < K ; i++) {
    for (int j = 0; j < d ; j++) {
      if (pca_dist(j) < abs(pca_data(i,j))) {
        pca_dist(j) = abs(pca_data(i,j));
      }
    }
  }

  for (int j = 0; j < d ; j++) {
    pca_volume = pca_volume + log(pca_dist(j));
  }
  //printf("LNC Call %d : %f, %f\n",index,pca_volume,rect_volume);
  correction = pca_volume - rect_volume;
  return correction;
}

//' mutual information function
//'
//' @param x numeric vector
//' @param low lower bound
//' @param high upper bound
//' @param res bounded numeric vector
//' @return bounded numeric vector
//' @export
// [[Rcpp::export]]
double rmi_mutual_information(Rcpp::NumericMatrix X,
                          Rcpp::NumericMatrix Y,
                          Rcpp::NumericMatrix XY,
                          int k,
                          int method,
                          Rcpp::NumericVector alpha) {

  int K     = k+1;
  int N     = X.nrow();
  int d_x   = X.ncol();
  int d_y   = Y.ncol();
  int d     = d_x + d_y;

  arma::mat  XY_dist(N,K);
  arma::imat XY_inds(N,K);


  arma::mat  X_dist(N,K);
  arma::imat X_inds(N,K);

  arma::mat  Y_dist(N,K);
  arma::imat Y_inds(N,K);

  double LNC_correction = 0;
  double proposed_correction = 0;

  if (method == 3 || method == 4) {
    get_nearest_neighbors(X, X_dist, X_inds, k);
    get_nearest_neighbors(Y, Y_dist, Y_inds, k);
  }

  double digamma_x = 0;
  double digamma_y = 0;
  double mi = boost::math::digamma(N) + boost::math::digamma(k);

  get_nearest_neighbors(XY, XY_dist, XY_inds, k);

  // //need these to be quere function
  // //get_nearest_neighbors(X, X_dist, X_inds, k);
  // //get_nearest_neighbors(Y, Y_dist, Y_inds, k);
  //
  // //Setting the flann index params
  // flann::IndexParams params;
  // params = flann::LinearIndexParams();
  //
  // // Setting flann search params
  // flann::SearchParams search_params;
  // //search_params.cores = 10;
  // //search_params.checks = 100;
  //
  // arma::mat X_query(d_x, N);
  // arma::mat X_temp_q(X.begin(), N, d_x, false);
  // X_query = arma::trans(X_temp_q); //transpose
  // flann::Matrix<double> X_ref_flann(X_query.memptr(), N, d_x);
  // flann::Index<flann::MaxDistance<double> > X_index(X_ref_flann, params);
  // X_index.buildIndex();
  //
  // arma::mat Y_query(d_y, N);
  // arma::mat Y_temp_q(Y.begin(), N, d_y, false);
  // Y_query = arma::trans(Y_temp_q); //transpose
  // flann::Matrix<double> Y_ref_flann(Y_query.memptr(), N, d_y);
  // flann::Index<flann::MaxDistance<double> > Y_index(Y_ref_flann, params);
  // Y_index.buildIndex();

  double epsilon   = 0;
  double epsilon_x = 0;
  double epsilon_y = 0;

  //printf("Finished Initialization\n");

  for (int i = 0; i < N; i++) {
   // printf("%d",i);
  // get epsilon x and espilon y
     epsilon = XY_dist(i,k-1);
   // printf("%f %d\n",epsilon,XY_inds(i,k-1));
    if (method == 1) {
      epsilon_x = epsilon;
      epsilon_y = epsilon;
      //printf("got here \n",i);
      int N_x = 0;
      int N_y = 0;

      for (int j = 0; j < N; j++) {
       // printf("%d\n",j);
        double dist_x = 0;
        double dist_y = 0;

        for (int m = 0; m < d_x; m++) {
          if (dist_x < abs(X(i,m)-X(j,m))) {
            dist_x = abs(X(i,m)-X(j,m));
          }
        }

        if (dist_x < epsilon_x) N_x++;

        for (int m = 0; m < d_y; m++) {
          if (dist_y < abs(Y(i,m)-Y(j,m))) {
            dist_y = abs(Y(i,m)-Y(j,m));
          }
        }
        if (dist_y < epsilon_y) N_y++;

      }
      //printf("eps: %f x: %f,%d y: %f,%d\n",epsilon, boost::math::digamma(N_x),N_x,boost::math::digamma(N_y),N_y);
      digamma_x = digamma_x + boost::math::digamma(N_x); //this KSG includes the point
      digamma_y = digamma_y + boost::math::digamma(N_y);

    } else { //KSG2

      epsilon_x = 0;
      epsilon_y = 0;

      for (int j = 1; j < K; j++) { //for each neighbor compute x and y dist

        double dist_x = 0;
        double dist_y = 0;

        for (int m = 0; m < d_x; m++) {
          if (dist_x < abs(X(i,m)-X(XY_inds(i,j),m))) {
            dist_x = abs(X(i,m)-X(XY_inds(i,j),m));
            //printf("dist_x %d: %f\n",i,dist_x);
          }
        }

        if (epsilon_x < dist_x) epsilon_x = dist_x;

        for (int m = 0; m < d_y; m++) {
          if (dist_y < abs(Y(i,m)-Y(XY_inds(i,j),m))) {
            dist_y = abs(Y(i,m)-Y(XY_inds(i,j),m));
            //printf("dist_y %d: %f\n",i,dist_y);
          }
        }

        if (epsilon_y < dist_y) epsilon_y = dist_y;
      }

      int N_x = 0;
      int N_y = 0;

      for (int j = 0; j < N; j++) {

        double dist_x = 0;
        double dist_y = 0;

        for (int m = 0; m < d_x; m++) {
          if (dist_x < abs(X(i,m)-X(j,m))) {
            dist_x = abs(X(i,m)-X(j,m));
          }
        }

        if (dist_x <= epsilon_x) N_x++;

        for (int m = 0; m < d_y; m++) {
          if (dist_y < abs(Y(i,m)-Y(j,m))) {
            dist_y = abs(Y(i,m)-Y(j,m));
          }
        }

        if (dist_y <= epsilon_y) N_y++;
        //printf("dx: %f, ex: %f, dy: %f, ey: %f\n",dist_x,epsilon_x,dist_y,epsilon_y);
      }

      //printf("x: %f,%f,%d y: %f,%f,%d\n",epsilon_x, boost::math::digamma(N_x),N_x,epsilon_y,boost::math::digamma(N_y),N_y);
      digamma_x = digamma_x + boost::math::digamma(N_x-1); //KSG2 does not include the point itself in the calculations
      digamma_y = digamma_y + boost::math::digamma(N_y-1);
    }

    if ( method == 3 || method == 4) {

      //joint space adjustment
      proposed_correction = compute_LNC(XY, XY_inds, i,d,K);

      if (proposed_correction < log(alpha(0))) {
       // printf("%d : %f, %f\n",i,proposed_correction,log(alpha(0)));
        LNC_correction = LNC_correction - proposed_correction;
      }

      if (method == 4) {
        if (d_x > 1) {
          proposed_correction = compute_LNC(X, X_inds, i,d_x,K);
          if (proposed_correction < alpha(1)) {
            LNC_correction = LNC_correction + proposed_correction;
          }
        }

        if (d_y > 1) {
          proposed_correction = compute_LNC(Y, Y_inds, i,d_y,K);
          if (proposed_correction < alpha(2)) {
            LNC_correction = LNC_correction + proposed_correction;
          }

        }

      }


    }

  }

  if (method == 1) {
    mi = mi - digamma_x/(double)N - digamma_y/(double)N;
  }
  if (method == 2) {
    mi = mi - digamma_x/(double)N - digamma_y/(double)N - 1/(double)k;
  }
  if (method == 3 || method == 4) {
    mi = mi - digamma_x/(double)N - digamma_y/(double)N - 1/(double)k + LNC_correction/(double)N;
  }
  return mi;
}

//' test get nearest function
//'
//' @param x numeric vector
//' @param low lower bound
//' @param high upper bound
//' @param res bounded numeric vector
//' @return bounded numeric vector
//' @export
// [[Rcpp::export]]
Rcpp::List test_get_nearest_neighbors(Rcpp::NumericMatrix data, int k) {
  arma::mat X_dist;
  arma::imat X_inds;
  get_nearest_neighbors(data, X_dist, X_inds, k);
  return Rcpp::List::create(Rcpp::Named("nn_dist") =  X_dist,Rcpp::Named("nn_inds") =  X_inds);
}





