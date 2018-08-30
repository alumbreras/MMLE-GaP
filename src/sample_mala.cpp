#include <RcppArmadillo.h>
#include <numeric>
#include "commons.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
SEXP sample_mala_cpp(const arma::ivec v_n, const arma::mat W, const arma::vec h_n_current,
                          double alpha = 1, double beta = 1, double delta = 0.01,
                          int iter = 100){
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif01(0.0,1.0);

  int F = W.n_rows;
  int K = W.n_cols;
  
  arma::mat eta_n_samples(K, iter);
  arma::vec eta_n_current(K);
  arma::vec eta_n(K);
  
  arma::vec grad_eta_current(K);
  arma::vec grad_eta(K);
  float a;
  
  eta_n_current = log(h_n_current);
  const arma::rowvec norms_W = sum(W, 0);
  
  
  grad_eta_current = grad_posterior_eta_cpp(eta_n_current, v_n, W, norms_W, alpha=alpha, beta=beta);
  
  double delta2 = pow(delta,2);
  for(int i=0; i<iter; i++){
    
    // Proposal
    for(int k=0; k<K; k++){
      std::normal_distribution<double> normal(eta_n_current[k] + 0.5 * delta2 * grad_eta_current[k], delta);
      eta_n[k] = normal(generator);
    }
    //eta_n = Rcpp::rnorm(K, eta_n_current + 0.5 * delta2 * grad_eta_current, delta);

    grad_eta = grad_posterior_eta_cpp(eta_n, v_n, W, norms_W, alpha=alpha, beta=beta);
    
    // Acceptance probability
    a = posterior_eta_cpp(eta_n, v_n, W, alpha, beta) -
      posterior_eta_cpp(eta_n_current, v_n, W, alpha, beta);

    double q_numer = - pow(norm(eta_n_current - eta_n - 0.5 * delta2 *  grad_eta), 2);
    double q_denom = - pow(norm(eta_n - eta_n_current - 0.5 * delta2 *  grad_eta_current), 2);
    a = a + 1/(2*delta2) * (q_numer - q_denom);
      
    // Accept or reject
    float random = unif01(generator);
    if(random < exp(a)){
      eta_n_current = eta_n;
      grad_eta_current = grad_eta;
    }
    
    // Store sample
    eta_n_samples.col(i) = eta_n_current;
  }
  // Transpose to have one sample per row
  eta_n_samples = eta_n_samples.t();
  
  // Back to the original space
  arma::mat h_samples = exp(eta_n_samples);
  
  return wrap(h_samples);
}
