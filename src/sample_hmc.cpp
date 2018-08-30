#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <numeric>
#include "commons.h"

using namespace Rcpp;
using namespace arma;

double U_cpp(arma::vec& eta_n, const arma::ivec& v_n, const arma::mat& W, 
             float alpha, float beta){
  return -posterior_eta_cpp(eta_n, v_n, W, alpha, beta);
}

arma::vec grad_U_cpp(const arma::vec& eta_n, const arma::ivec& v_n, 
                     const arma::mat& W, const arma::rowvec& norms_W,
                     float alpha, float beta){
  
  return - grad_posterior_eta_cpp(eta_n, v_n, W, norms_W, alpha, beta);
}

// Samples using the logtransform
// [[Rcpp::export]]
SEXP sample_hmc_cpp(const arma::ivec& v_n, const arma::mat& W, arma::vec h_n_current,
                    float alpha = 1, float beta = 1,
                    int L=10, float epsilon = 0.01,
                    int iter=100){
  
  int K = W.n_cols;
  int F = W.n_rows;

  // Pre-compute column norms
  arma::irowvec u(F);
  u.ones();
  arma::rowvec norms_W = u * W;

  arma::mat h_n_samples(K, iter);   // traces of p
  arma::vec q(K);                   // position
  arma::vec current_q(K);
  arma::vec p(K);                   // momentum
  arma::vec current_p(K);

  current_q = log(h_n_current);
  h_n_samples.col(0) = current_q;
  for(int i=1; i<iter; i++){
    q = current_q;
    p = rnorm(K, 0,1); // independent standard normal variates
    current_p = p;

    // Make a half step for momentum at the beginning
    p = p - epsilon * grad_U_cpp(q, v_n, W, norms_W, alpha, beta) / 2;

    // Alternate a full step for the position
    for(int l=0; l<L; l++){
      // Make a full step for the position
      q = q + epsilon * p;
      //q = abs(q); // Bouncing. Recommended by Nico.

      // Make a full step for the momentum, except at the end of trajectory
      if(l != L) {
        p = p - epsilon * grad_U_cpp(q, v_n, W, norms_W, alpha, beta);
      }
    }
  
    // Make a half step for momentum at the end.
    p = p - epsilon * grad_U_cpp(q, v_n, W, norms_W, alpha=alpha, beta=beta) / 2;

    // Negate momentum at end of trajectory to make the proposal symmetric
    p = -p;
    
    //Evaluate potential and kinetic energies at start and end of trajectory
    double current_U  = U_cpp(current_q, v_n, W, alpha, beta);
    double current_K = as_scalar(current_p.t() * current_p) / 2 ;
    double proposed_U = U_cpp(q, v_n, W, alpha, beta);
    double proposed_K = as_scalar(p.t() * p) / 2;
    
    //Accept or reject the state at the end of trajectory, returning either
    // the position at the end of the trajectory or the initial position
    double current_H = current_U + current_K;
    double proposed_H = proposed_U  + proposed_K;
    double uni = ::Rf_runif (0, 1);
    double accept = exp(current_H -  proposed_H);
    if(uni < accept){
      current_q = q;  // accept
    }

    h_n_samples.col(i) = current_q;
  }
  
  h_n_samples = exp(h_n_samples.t());
  return(wrap(h_n_samples));
}