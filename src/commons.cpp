// Functions that may be used by different .cpp files
// Thanks to the export instruction, they can also be used from the R code
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double posterior_h_cpp(const arma::vec& h_n, const arma::ivec& v_n, const arma::mat& W, 
                        float alpha, float beta){
  double logp = 0;
  
  if(any(h_n < 0)){
    return - arma::datum::inf;
  }
  
  int K = W.n_cols;
  int F = W.n_rows;
  
  arma::dvec lambdas = W * h_n;
  
  for(int f = 0; f < F; f++){
    logp += v_n[f] * log(lambdas[f]) - lambdas[f]; //- std::lgamma(v_n[f]+1);
  }

  for(int k = 0; k < K; k++){
    logp += alpha * log(h_n[k]) - log(h_n[k]) - beta*h_n[k];
  }
  
  return logp;
}


double posterior_h_cpp_3(const arma::vec h_n, const arma::ivec& v_n, const arma::mat& W, 
                      float alpha, float beta){
  
  double logp = 0;
  
  if(any(h_n < 0)){
    return - arma::datum::inf;
  }
  
  int K = W.n_cols;
  int F = W.n_rows;
  
  arma::dvec lambdas = W * h_n;
  logp = sum(v_n % log(lambdas)) - sum(lambdas) + 
          sum((alpha-1)*log(h_n) - beta*h_n);
  
  return logp;
}


// [[Rcpp::export]]
double posterior_eta_cpp(const arma::vec& eta_n, const arma::ivec& v_n, const arma::mat& W, 
                        float alpha, float beta){
  arma::vec exp_eta = exp(eta_n);
  double logp = posterior_h_cpp(exp_eta, v_n, W, alpha, beta) + sum(eta_n);
  return logp;
}


arma::vec grad_posterior_eta_cpp(const arma::vec& eta_n, const arma::ivec& v_n, 
                                 const arma::mat& W, const arma::rowvec& norms_W,
                                 float alpha, float beta){
  
  int K = eta_n.size();  
  arma::vec dh = zeros(K);
  arma::vec lambdas = W * exp(eta_n);
  for(int k=0; k < K; k++){
    dh[k] = alpha - exp(eta_n[k])*(beta + norms_W[k]) + sum(v_n % (W.col(k)*exp(eta_n[k])/lambdas));
  }
  return dh;
}

// [[Rcpp::export]]
SEXP grad_posterior_eta_cpp_R(const arma::vec& eta_n, const arma::ivec& v_n, 
                                 const arma::mat& W, const arma::rowvec& norms_W,
                                 float alpha, float beta){
  
  arma::vec dh = grad_posterior_eta_cpp(eta_n, v_n, W, norms_W, alpha, beta);
  return(wrap(dh));
}



// Multivariate gaussian samples
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

//Exact computation of p(V|W)
//double likelihood(arma:vec v, arma:mat W){
//}
