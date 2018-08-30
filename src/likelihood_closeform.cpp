#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <numeric>

using namespace Rcpp;
using namespace arma;

double likelihood_C_W(const	arma::imat& C, 
                      const arma::mat& W,
                      float alpha, float beta){

  int K = W.n_cols;
  int F = W.n_rows;
  
  arma::colvec denominators = beta + sum(W,0).t();
  arma::ivec C_norms = sum(C,0).t();
  
  double prob = 0;
  for (int k=0; k<K; k++){
    double prob_ = 0;
    
    prob_ += std::lgamma(alpha + C_norms[k]) - std::lgamma(alpha);
    for(int f=0; f<F; f++){
      prob_ -= std::lgamma(C(f,k)+1);
    }
    prob = prob + prob_;
    prob += log(pow((beta/denominators[k]), alpha));

    prob_ = 0;
    for(int f=0; f<F; f++){
      prob_ += log(pow((W(f,k)/denominators[k]), C(f,k)));	  		
    }

    prob = prob + prob_;
  }
  return prob;
}

double likelihood_recursive(const arma::icolvec& v, 
                            arma::imat& C, 
                            const arma::mat& W, 
                            float alpha, float beta,
                            int f, int pos, int remaining) {
  double like = 0;
  
  // If there is no remaining left, we completed a new composition for v[f]
  if (remaining == 0) { 
    
    // If features left, get the combinations of v[f+1]
    if(f < (C.n_rows-1)){
      like = likelihood_recursive(v, C, W, alpha, beta, f+1, 0, v[f+1]);
      return like;
    } 
    // If last f, then we are done and we completed a new C
    else {
      // compute likelihood
      like = std::exp(likelihood_C_W(C, W, alpha, beta)); // Ned non log because we have to \sum_C p(C|W)
      return like;
    }
  }
  
  // If position pointer got out of the vector, 
  // then there is nothing to do
  if (pos == C.n_cols) { return 0; }
  
  // Else, continue allocating the remaining in all possible ways
  for (int i = remaining; i >= 0; --i) {
    C(f, pos) = i;
    like += likelihood_recursive(v, C, W, alpha, beta, f, pos + 1, remaining - i);
  }
  return like;
}

// [[Rcpp::export]]
double likelihood_V_W_cpp(const arma::imat& V,  const arma::mat& W, float alpha, float beta) {
  
  int N = V.n_cols;
  int F = V.n_rows;
  int K = W.n_cols;
  arma::imat C(F,K);
  C.zeros();
  
  // Compute likelihood column by column
  double like = 0;
  for(int n=0; n<N; n++){
    like += std::log(likelihood_recursive(V.col(n), C, W, alpha, beta, 0, 0, V(0,n)));
  }
  
  return like;
}