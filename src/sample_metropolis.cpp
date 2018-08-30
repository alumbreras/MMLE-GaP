#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <numeric>
#include "commons.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP sample_metropolis_h_reparam_cpp(arma::ivec v_n, const arma::mat& W, arma::vec h_n_current,
                      float alpha = 1, float beta = 1, float step = 1, int iter=100){

  int F = W.n_rows;
  int K = W.n_cols;
  
  arma::mat eta_n_samples(K, iter);
  arma::vec eta_n_current = log(h_n_current);
  arma::vec eta_n;
  float a;
  
  eta_n_samples.col(0) = eta_n_current;
  for(int i=1; i<iter; i++){

    // Proposal
    arma:vec jump = Rcpp::rnorm(K, 0, step);
    eta_n = eta_n_current + jump;
  
    // Acceptance probability
    a = posterior_eta_cpp(eta_n, v_n, W, alpha, beta) -
        posterior_eta_cpp(eta_n_current, v_n, W, alpha, beta);
    
    // Accept or reject
    if(R::runif(0, 1) < exp(a)){
      eta_n_current = eta_n;
    }
    
    // Store sample
    eta_n_samples.col(i) = eta_n_current;
  }
  
  // Transpose to have one sample per row
  eta_n_samples = eta_n_samples.t();
  
  // Back to the original space
  arma::mat h_samples = exp(eta_n_samples);

  return(wrap(h_samples));
}

