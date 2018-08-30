// Samples p(C | V, W) without using H
// It considers that c_{fkn} is the number of counts f in topic k and document n,
// and resamples the topic assigment of each count that contributes to V.
// (see, for instance, Buntine's paper "Discrete Component Analysis")

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <numeric>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
SEXP sample_gibbs_z_cpp(const arma::vec& v_n, 
                        const arma::mat& W, 
                        arma::imat C,
                        double alpha = 1, double beta = 1, int iter=100,
                        double burnin = 0.5){
  
  int K = W.n_cols;
  int F = W.n_rows;
  int nburnin = std::floor(iter*burnin);
  int nsamples = iter - std::floor(iter*burnin);
  
  // Rcout << "=============sample_gibbs_z_cpp============" << std::endl;
  // Rcout << "v_n" << std::endl;
  // Rf_PrintValue(wrap(v_n));
  // Rcout << "C" << std::endl;
  // Rf_PrintValue(wrap(C));
  // Rcout << "W" << std::endl;
  // Rf_PrintValue(wrap(W));
  
  arma::mat    C_samples_mean(F, K);
  C_samples_mean.zeros();
  int j = 0;
  
  arma::irowvec     L(K);
  arma::rowvec      probs(K);
  arma::imat        C_initial(F,K);
  arma::irowvec     ans(K);
  
  arma::rowvec denominators;
  denominators = beta + sum(W, 0);
  
  L = arma::sum(C, 0); // tells how many counts in each k
  
  C_samples_mean = C_samples_mean + C; //init sample
  j = j + 1;
  
  // Resample the topic assigment k of each count
  for(int i=1; i<iter; i++){
    
    // Fix the original matrix so that we can sample each occurrence
    for(int ii=0; ii<F; ii++){
      for(int jj=0; jj<K; jj++){
        C_initial(ii,jj) = C(ii,jj);
      }
    }
    // Rcout << "C_initial"  << std::endl;
    // Rf_PrintValue(wrap(C_initial));
    for(int f=0; f<F; f++){
      
      // v_n[f] total counts of feature f
      // for each occurent of f in the document and topic k, re-sample its topic
      for(int k=0; k<K; k++){
          if(C_initial(f,k) == 0) { 
            continue;
          }
         // Use a C_initial snapshot so that each occurrence is re-sampled only once
         for(int l=0; l < C_initial(f,k); l++){
             L[k] = L[k] - 1;    // remove occurrence
             C(f, k) = C(f, k) - 1; // remove occurrence
        
           // Sample new k for the occurrence
            probs = (alpha + L) % (W.row(f)/denominators);
            probs = probs/sum(probs);
            ans = ans.zeros();
            rmultinom(1, probs.begin(), K, ans.begin());
        
            C.row(f) = C.row(f) + ans;  // add occurrence
            L = L + ans;     // add occurrence
        }
      }
    }
    // After burning, store samples
    // (or sufficient statistics, if samples won't be needed)
    if(i >= nburnin){
      C_samples_mean = C_samples_mean + C; 
      j = j + 1;
    }
    C_samples_mean = C_samples_mean/j;
  }
  
  return Rcpp::wrap(C_samples_mean);
}