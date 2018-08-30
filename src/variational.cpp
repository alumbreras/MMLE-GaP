#include <RcppArmadillo.h>
#include <numeric>
#include <boost/math/special_functions/digamma.hpp>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
SEXP variational_ch_cpp(arma::ivec v_n, arma::mat W, 
                        arma::vec alpha_var, arma::vec beta_var, 
                        double alpha = 1, double beta = 1, 
                        int maxiters = 100){
  
  // Containers
  int K = W.n_cols;
  int F = W.n_rows;
  arma::mat prob_var(F,K);
  arma::mat E_c(F,K);
  arma::vec E_ln_h(K);  
  arma::dvec lp(maxiters);
  
  // Precomputations
  arma::irowvec rowones(K);
  rowones.ones();
  arma::imat v_n_template = v_n * rowones;  
  arma::vec W_norms = sum(W, 0).t();
  
  
  
  for(int i=0; i<maxiters; i++){
    //Rcpp::Rcout << "iteration: " <<  i << std::endl;
    
    // Compute sufficient statistics
    for(int k=0; k<K; k++){
      E_ln_h[k] = boost::math::digamma(alpha_var[k]) - log(beta_var[k]);
    }

    // Update variational parameters for q(C)
    arma::vec Wh = W * exp(E_ln_h);
    for(int f=0; f<F; f++){
      for(int k=0; k<K; k++){
        prob_var(f,k) = W(f,k) * exp(E_ln_h[k])/Wh[f];
      }
    }
     
    // Compute sufficient statistics
    E_c = prob_var % v_n_template;
     
    // Update variational parameters for q(H)
    arma::vec E_c_norms = sum(E_c, 0).t();
    alpha_var = alpha + E_c_norms; 
    beta_var  = beta + W_norms;
    arma::vec E_h(K);
    E_h = alpha_var/beta_var; 
    
    // Compute lower bound
    double lgammma_alpha = std::lgamma(alpha) * K;
    arma::vec lgamma_alpha_var(K);
    for(int k=0; k<K; k++){
      lgamma_alpha_var[k] = std::lgamma(alpha_var[k]);
    }
    
    double sum_logfac_v = 0;
    for(int f=0; f<F; f++){
      sum_logfac_v += std::lgamma(v_n[f]+1);

      
    }
    
    double epsilon = 0.00001; // to avoid 0*-inf

    double lp1 = - accu(W * E_h) + accu(E_c % (log(W+epsilon) + repmat(E_ln_h.t(), F, 1)));
    //Rcout << "lp1: " <<  lp1 << std::endl;
    
    double lp2 = sum((alpha-1)*E_ln_h - beta * E_h + alpha * log(beta)) - lgammma_alpha;
    //Rcout << "lp2: " <<  lp2 << std::endl;
    
    double lp3 = sum_logfac_v + accu(E_c % log(prob_var+epsilon));
    
    double lp4 = sum(alpha_var % log(beta_var) - lgamma_alpha_var + (alpha_var-1) % E_ln_h - alpha_var);
    //Rcout << "lp4: " <<  lp4 << std::endl;
    
    lp[i] = lp1 + lp2 - lp3 - lp4;
    //std::cout << "lp: " <<  lp[i] << std::endl;
    
  }
    
  return Rcpp::List::create(Rcpp::Named("q_h_alpha") = wrap(alpha_var),
                            Rcpp::Named("q_h_beta")  = wrap(beta_var),
                            Rcpp::Named("lp") = wrap(lp));
  
}