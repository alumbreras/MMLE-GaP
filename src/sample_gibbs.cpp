#include "configs.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <numeric>

using namespace Rcpp;
using namespace arma;

arma::irowvec sample_gibbs_c_cpp(const int v_fn, 
                                 const arma::rowvec& w_f, 
                                 const arma::vec& h_n){

  int K = h_n.size();
  arma::irowvec c = arma::zeros<irowvec>(K);
  
  // If v_fn is zero, then its corresponding c is all zero's.
  if(v_fn == 0){
    return c;
  }
  
  // Compute probabilities of c_1,...c_k
  arma::vec probs(K);
  for(int k=0; k<K; k++){
    probs[k] = w_f[k] * h_n[k];  
  }
  
  double normalization = arma::accu(probs);
  //if(normalization < 0.0001){
  //  return c;
  //}
  probs = probs/normalization;

  // Draw v_fn experiments to fill the vector c_1,...c_k
  rmultinom(v_fn, probs.begin(), K, c.begin());
  return(c);
}

double sample_gibbs_h_cpp(const arma::vec& W_k, const arma::ivec& C_k, 
                          const double alpha, const double beta){
  double alpha_post = alpha + arma::sum(C_k); //std::accumulate(C_k.begin(), C_k.end(), 0.0);
  double beta_post = 1/(beta + arma::sum(W_k)); //std::accumulate(W_k.begin(), W_k.end(), 0.0));
  return R::rgamma(alpha_post, beta_post);
}


// [[Rcpp::export]]
SEXP sample_gibbs_cpp(const arma::vec& v_n, const arma::mat& W, arma::vec h_n,
                      double alpha = 1, double beta = 1, int iter=100,
                      double burnin = 0.5){

  int K = W.n_cols;
  int F = W.n_rows;
  int nburnin = std::floor(iter*burnin);
  int nsamples = iter - std::floor(iter*burnin);
  int j = 0;
  
  // We do not need the iters dimensions in the C_samples.
  // we can just sum over iterations, and this is all we need 
  // for the M-step
  arma::mat    C_samples_mean(F, K);
  arma::vec    h_samples_mean(K);
  C_samples_mean.zeros();
  h_samples_mean.zeros();
  
  arma::imat   C(F,K);
  arma::mat    h_samples(K, nsamples);
  
  // First sample with h_init
  for(int f=0; f<F; f++){
   C.row(f) = sample_gibbs_c_cpp(v_n[f], W.row(f), h_n);
  }

  //If no burning, the initial sample is stored
  if(nburnin == 0){
    C_samples_mean = C_samples_mean + C;
    h_samples_mean = h_samples_mean + h_n;
    h_samples.col(0) = h_n;
    j = j + 1;
  }
  
  for(int i=1; i<iter; i++){
    for(int f=0; f<F; f++){
       C.row(f) = sample_gibbs_c_cpp(v_n[f], W.row(f), h_n);
    }

    for(int k=0; k<K; k++){
      h_n[k] = sample_gibbs_h_cpp(W.col(k), C.col(k), alpha, beta);
    }
    
    // After burning, store samples
    // (or sufficient statistics, if samples won't be needed)
    if(i >= nburnin){
      C_samples_mean = C_samples_mean + C;
      h_samples_mean = h_samples_mean + h_n;
      h_samples.col(j) = h_n;
      j = j + 1;
    }
  }
  C_samples_mean = C_samples_mean/j;
  h_samples_mean = h_samples_mean/j;
  
  h_samples = h_samples.t();

  // We might try to permute also C dimensions to return (iter, F,K)
  // But this is not trivial. Do it with R functions, even from here
   return Rcpp::List::create(Rcpp::Named("C_n_samples_mean") = wrap(C_samples_mean),
                             Rcpp::Named("h_n_samples_mean") = wrap(h_samples_mean),
                             Rcpp::Named("h_n_samples") = wrap(h_samples));
}

// The same Gibbs sampler, but looping through all N 
// so that we do not need to do it in R
// [[Rcpp::export]]
SEXP sample_gibbs_cpp_N(const arma::sp_mat& V, const arma::mat& W, arma::mat H_init,
                       double alpha = 1, double beta = 1, int iter=100,
                       double burnin = 0.5){

  int K = W.n_cols;
  int F = W.n_rows;
  int N = V.n_cols;
  int nburnin = std::floor(iter*burnin);
  int nsamples = iter - std::floor(iter*burnin);
  
  //arma::ivec    v_n(F);
  
  //Objects to be returned to R
  arma::cube   C_samples_mean_N(F, K, N);
  arma::mat    H_samples_mean_N(K, N);
  arma::cube   H_samples_N(nsamples, K, N);
  C_samples_mean_N.zeros();
  H_samples_mean_N.zeros();
  
  //Objects to be re-filled at each n or iteration
  arma::vec    h_n(K);
  arma::imat   C(F,K);
  arma::mat    h_samples(K, nsamples);
  arma::mat    C_samples_mean(F, K);
  arma::vec    h_samples_mean(K);  
  
  int j = 0;
  
  for(int n=0; n<N ; n++){
    if(!((n+1) % 1000)){
      Rcout << "\n gibbs in column: " << n << std::endl;
    }
    
    j = 0;
  
    C_samples_mean.zeros();
    h_samples_mean.zeros();  
    h_samples.zeros();
    
    h_n = H_init.col(n);
    //v_n = V.col(n);
    
    // First sample with h_init
    for(int f=0; f<F; f++){
      C.row(f) = sample_gibbs_c_cpp(V(f,n), W.row(f), h_n);
    }
    
    //If no burning, the initial sample is stored
    if(nburnin == 0){
      C_samples_mean = C_samples_mean + C;
      h_samples_mean = h_samples_mean + h_n;
      h_samples.col(0) = h_n;
      j = j + 1;
    }
    
    for(int i=1; i<iter; i++){    
  
      // First sample with h_init
      for(int f=0; f<F; f++){
        C.row(f) = sample_gibbs_c_cpp(V(f,n), W.row(f), h_n);
      }
      for(int k=0; k<K; k++){
        h_n[k] = sample_gibbs_h_cpp(W.col(k), C.col(k), alpha, beta);
      }
      
      // After burning, store samples
      // (or sufficient statistics, if samples won't be needed)
      if(i >= nburnin){
        C_samples_mean = C_samples_mean + C;
        h_samples_mean = h_samples_mean + h_n;
        h_samples.col(j) = h_n;
        j = j + 1;
      }
    }
    
    C_samples_mean = C_samples_mean/j;
    h_samples_mean = h_samples_mean/j;
    
    H_samples_N.slice(n) =  h_samples.t();
    C_samples_mean_N.slice(n) = C_samples_mean;
    H_samples_mean_N.col(n) = h_samples_mean;
  }

  return Rcpp::List::create(Rcpp::Named("C_samples_mean") = wrap(C_samples_mean_N),
                            Rcpp::Named("H_samples_mean") = wrap(H_samples_mean_N),
                            Rcpp::Named("H_samples") = wrap(H_samples_N));    
}
  