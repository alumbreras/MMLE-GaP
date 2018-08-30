#include <RcppArmadillo.h>
#include <numeric>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <iostream>
//#include "commons.h"

namespace S {

  // Stores position q and momentum p
  struct pq_struct {
    arma::vec q;
    arma::vec p;
  };
  
  // Stores position-momentum of forward and backwards path
  struct Tree {
    arma::vec minus_q;
    arma::vec plus_q;
    arma::vec minus_p;
    arma::vec plus_p;
    arma::vec q_prima;
    int n;
    int s;
  };

}

using namespace Rcpp;
using namespace arma;

// float posterior_h_cpp(arma::vec& h_n, const arma::ivec& v_n, const arma::mat& W, 
//                       float alpha, float beta){
//   
//   float logp = 0;
//   
//   if(any(h_n < 0)){
//     return - arma::datum::inf;
//   }
//   
//   int K = W.n_cols;
//   int F = W.n_rows;
//   
//   arma::vec lambdas = W * h_n;
//   
//   for(int f = 0; f < F; f++){
//     logp += v_n[f] * log(lambdas[f]) - lambdas[f]; //- log(std::tgamma(v_n[f]+1));  
//   }
//   for(int k = 0; k < K; k++){
//     logp += (alpha-1)*log(h_n[k]) - beta*h_n[k];
//   }
//   
//   // Similar alternative:
//   //logp = sum(v_n % log(lambdas)) - beta*sum(lambdas) + 
//   //        sum((alpha-1)*log(h_n) - beta*h_n);
//   
//   return logp;
// }

// [[Rcpp::export]]
double posterior_h_cpp_(arma::vec& h_n, const arma::ivec& v_n, const arma::mat& W, 
                      float alpha, float beta){
  // uses pow() so that 0^0=0 rather than 0*log(0) = NaN
  // this mimics a Gamma*Poisson computation in R
  
  double logp = 0;
  
  if(any(h_n < 0)){
    return - arma::datum::inf;
  }
  
  int K = W.n_cols;
  int F = W.n_rows;
  
  arma::dvec lambdas = W * h_n;
  
  for(int f = 0; f < F; f++){
    //logp += v_n[f] * log(lambdas[f]) - lambdas[f]; //- log(std::tgamma(v_n[f]+1));
    logp += log(pow(lambdas[f], v_n[f])) - lambdas[f];
  }

  for(int k = 0; k < K; k++){
    //logp += (alpha-1)*log(h_n[k]) - beta*h_n[k];
    logp += log(pow(h_n[k],(alpha-1))) - beta*h_n[k];
  }
  
  // Similar alternative:
  //logp = sum(v_n % log(lambdas)) - beta*sum(lambdas) + 
  //        sum((alpha-1)*log(h_n) - beta*h_n);
  
  return logp;
}

// Posterior wrapper
double loglike_cpp(arma::vec& h_n, const arma::ivec& v_n, const arma::mat& W, 
                   float alpha, float beta){
  return posterior_h_cpp_(h_n, v_n, W, alpha, beta);
}

// Gradient of the log posterior
arma::vec grad_loglike_cpp(arma::vec h_n, arma::ivec v_n, 
                           const arma::mat& W, const arma::rowvec& norms_W,
                           float alpha, float beta){
  
  int K = h_n.size();  
  arma::vec dh = zeros(K);
  arma::vec lambdas = W * h_n; // this is the slowest bottleneck
  for(int k=0; k < K; k++){
    dh[k] = (alpha-1)/h_n[k] - (beta + norms_W[k]) + sum(v_n % (W.col(k)/lambdas));
  }
  return dh;
}


// arma::vec grad_loglike_cpp4(arma::vec h_n, arma::ivec v_n, 
//                       const arma::mat& W, const arma::vec& norms_W,
//                       float alpha, float beta){
//   
//   int K = h_n.size();
//   int F = W.n_rows;
//   
//   arma::vec dh = zeros(K);
//   arma::vec lambdas = W * h_n;
//   for(int k=0; k < K; k++){
//     dh[k] = (alpha-1)/h_n[k] - (beta + norms_W[k]);
//     for(int f=0; f < F; f++){
//       dh[k] += v_n[f] * W(f,k) / lambdas[f];
//     }
//   }
//   return dh;
// }


// Performs one leapfrom step (NUTS paper, Algorithm 1)
S::pq_struct leapfrog(arma::vec q, arma::vec p, float epsilon, 
                   const arma::ivec v_n, const arma::mat& W, const arma::rowvec& norms_W, 
                   float alpha, float beta){
  
  S::pq_struct pq;
  
  p += epsilon * 0.5 * grad_loglike_cpp(q, v_n, W, norms_W, alpha=alpha, beta=beta);
  q += epsilon * p;
  q  = abs(q); // Bouncing. Recommended by Nico (to stay in non-negative values)
  p += epsilon * 0.5 * grad_loglike_cpp(q, v_n, W, norms_W, alpha=alpha, beta=beta);
  
  pq.q = q; pq.p = p;
  return pq;
}

S::Tree BuildTree(arma::vec q, arma::vec p, float u, int v, int j, float epsilon,
               arma::ivec v_n, const arma::mat& W, const arma::rowvec& norms_W, 
               float alpha, float beta){
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif01(0.0,1.0);
  
  S::Tree tree;
  arma::vec q_prima;
  arma::vec p_prima;
  arma::vec minus_q;
  arma::vec minus_p;
  arma::vec plus_q;
  arma::vec plus_p;
  int n_prima;
  int s_prima;
  arma::vec q_prima2;
  int n_prima2;
  int s_prima2;
  
  float delta_max = 0;
  
  if(j==0){
    // Base case - take one leapfrog step in each direction
    // Metropolis accept-reject will be done in the main function
    S::pq_struct pq = leapfrog(q, p, v*epsilon, v_n, W, norms_W, alpha, beta);
    q_prima = pq.q;
    p_prima = pq.p;
    float U_prima = loglike_cpp(q_prima, v_n, W, alpha, beta); 
    float K_prima = 0.5 * as_scalar(p_prima.t() * p_prima);
    int n_prima   = u <= exp(U_prima - K_prima);
    int s_prima   = (U_prima - K_prima ) > (log(u)-delta_max);
    
    tree.minus_q = q_prima;
    tree.minus_p = p_prima;
    tree.plus_q  = q_prima;
    tree.plus_p  = p_prima;
    tree.q_prima = q_prima;
    tree.n = n_prima;
    tree.s = s_prima;
    return tree;
    
  } else {
    // Recursion -- implicitly build the left and right subtrees
    S::Tree tree = BuildTree(q, p, u, v, j-1, epsilon, v_n, W, norms_W, alpha, beta);
    minus_q = tree.minus_q;
    minus_p = tree.minus_p;
    plus_q  = tree.plus_q;
    plus_p  = tree.plus_p;
    q_prima = tree.q_prima;
    n_prima = tree.n;
    s_prima = tree.s;
    
    if(s_prima == 1){
      if(v == -1){
        S::Tree tree = BuildTree(minus_q, minus_p, u, v, j-1, epsilon, v_n, W, norms_W, alpha, beta);
        minus_q  = tree.minus_q;
        minus_p  = tree.minus_p;
        q_prima2 = tree.q_prima;
        n_prima2 = tree.n;
        s_prima2 = tree.s;
      } else {
        S::Tree tree  = BuildTree(plus_q, plus_p, u, v, j-1, epsilon, v_n, W, norms_W, alpha, beta);
        plus_q   = tree.plus_q;
        plus_p   = tree.plus_p;
        q_prima2 = tree.q_prima;
        n_prima2 = tree.n;
        s_prima2 = tree.s;
      }
      
      if(unif01(generator) < n_prima2/(n_prima + n_prima2)){
        q_prima = q_prima2; 
      }
      
      arma::rowvec diff_q = (plus_q - minus_q).t();
      int uturn = (as_scalar(diff_q * minus_p) >= 0) && (as_scalar(diff_q * plus_p) >= 0);
      s_prima  =  s_prima2 * uturn;
      n_prima +=  n_prima2;
    }
    
    tree.minus_q = minus_q;
    tree.minus_p = minus_p;
    tree.plus_q  = plus_q;
    tree.plus_p  = plus_p;
    tree.q_prima = q_prima;
    tree.n = n_prima;
    tree.s = s_prima;
    return tree;
  }
}

// [[Rcpp::export]]
SEXP sample_nuts_cpp(arma::ivec v_n, const arma::mat& W, arma::vec h_n_current,
                          double alpha = 1, double beta = 1,
                          float epsilon = 0.01,
                          int iter=100){
  
  int K = W.n_cols;
  int F = W.n_rows;
  
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif01(0.0,1.0);
  
  // Pre-compute column norms
  arma::rowvec u_(F);
  arma::rowvec norms_W = u_ * W;
  
  //Rf_PrintValue(wrap(norms_W));
  arma::mat h_n_samples(K, iter);  // traces of p
  arma::vec q(K);                 // position
  arma::vec current_q(K);
  arma::vec p(K);                 // momentum
  arma::vec current_p(K);
  arma::vec minus_q(K);
  arma::vec minus_p(K);
  arma::vec plus_q(K);
  arma::vec plus_p(K);
  arma::vec q_prima(K);
  int n_prima;
  int s_prima;
  current_q = h_n_current;
  h_n_samples.col(1) = h_n_current;
  
  //cout << "Start iters " << std::endl;
  for(int i=1; i<iter; i++){
    Rcout << "Sample: " << i << std::endl;
    
    // Sample new momentum (K independent standard normal variates)
    std::normal_distribution<double> normal(0,1);
    for(int k=1; k<K; k++){
      p[k] = normal(generator);
    }
    
    // Draw slice variable  
    float limit_sup = exp( loglike_cpp(current_q, v_n, W, alpha, beta) - 0.5* as_scalar(p.t()*p));
    std::uniform_real_distribution<double> distribution(0.0,limit_sup);
    float u = distribution(generator);
    
    // Initialize position and momentum.
    minus_q = current_q;  // position in the backward path
    plus_q = current_q;   // position in the forward path
    minus_p = p;          // momentum in the backward path
    plus_p = p;           // momentum in the forward path
    current_p = p;        // position from where backward and forward paths start
    int j = 0;            // Step number
    int n = 1;
    int s = 1;            // Binary. 1 while no U-turn detected
    int v;                // Coin to decide direction
    
    // While no U-turn
    while(s==1){
      
      // Draw direction
      if(unif01(generator) < 0.5) { v = -1;} else { v = 1;}
      
      
      if(v == -1){ // if backwards
        S::Tree tree = BuildTree(minus_q, minus_p, u, v, j, epsilon, v_n, W, norms_W, 
                              alpha, beta);
        minus_q = tree.minus_q;
        minus_p = tree.minus_p;
        q_prima = tree.q_prima;
        n_prima = tree.n;
        s_prima = tree.s;
      } else {       // if forwards
        S::Tree tree = BuildTree(plus_q,   plus_p, u, v, j, epsilon, v_n, W, norms_W, 
                              alpha, beta);
        plus_q  = tree.plus_q;
        plus_p  = tree.plus_p;
        q_prima = tree.q_prima;
        n_prima = tree.n;
        s_prima = tree.s;
      }
      
      if(s_prima == 1){ 
        // Accept - Reject
        if(unif01(generator) < std::min(1, n_prima/n)){ 
          current_q = q_prima; // Accept proposal
        }
      }
      
      n += n_prima;
      arma::rowvec diff_q = (plus_q - minus_q).t();
      s = s_prima * (as_scalar(diff_q * minus_p) >= 0) && (as_scalar(diff_q * plus_p) >= 0);
      j++; 
      Rcout << "Depth j " << j << std::endl;
    } // end while
    
    h_n_samples.col(i) = current_q;
    
  } // end for
  h_n_samples = h_n_samples.t();
  return(Rcpp::wrap(h_n_samples));
  
}