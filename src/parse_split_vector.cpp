#include <RcppArmadillo.h>
using namespace std;

// [[Rcpp::export]]
void parse_split_vector(Rcpp::NumericVector splits,
                        arma::icolvec& d,
                        arma::icolvec& start_d,
                        arma::icolvec& end_d) {
  int vars = d.n_elem;
  // this code is nasty, I wished it was clearer
  start_d[0] = 0; // c(1,2,3) -> 0,0 - 1,2 - 3,5 - 0,1,2,3,4,5
  start_d[1] = 0;
  end_d[1]   = splits[0]-1;
  for (int i = 1; i < vars; i++) {
    d[i]        =  splits[i-1];
    start_d[i+1] = start_d[i] + d[i];
    end_d[i]   = start_d[i+1]  - 1;
    d[0]       = d[0] + d[i];
  }
  end_d[0] = end_d[vars-1];
}
