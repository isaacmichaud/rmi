#include <RcppArmadillo.h>
using namespace std;

double lnc_compute(const arma::mat  & data,
                   const arma::imat & nn_inds,
                   int index,
                   int start_ind,
                   int end_ind) {
  // initial declarations
  int K = nn_inds.n_cols;
  int d = (end_ind - start_ind) + 1;
  arma::colvec mean = arma::zeros(d);
  arma::mat X       = arma::zeros(K,d);
  double correction = 0;

  // center the data
  for (int i = 0; i < K ; i++) {
    for (int j = 0; j < d ; j++) {
      mean(j) = mean(j) + data(nn_inds(index,i),j)/K;
    }
  }

  for (int i = 0; i < K ; i++) {
    for (int j = 0; j < d ; j++) {
      X(i,j) = data(nn_inds(index,i),j) - mean(j);
    }
  }

  // compute coordinate aligned volume
  double rect_volume = 0;
  arma::colvec rect_dist = arma::zeros(d);
  for (int i = 0; i < K ; i++) {
    for (int j = 0; j < d ; j++) {
      if (rect_dist(j) < abs(X(i,j))) {
        rect_dist(j) = abs(X(i,j));
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
  arma::svd(U, s, V, X);

  arma::mat pca_data(K,d);
  pca_data = X*V;

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








