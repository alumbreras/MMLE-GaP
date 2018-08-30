#ifndef COMMONS_INCLUDED
#define COMMONS_INCLUDED

// Functions that may be used by different .cpp files
double posterior_h_cpp(const arma::vec& h_n, const arma::ivec& v_n, const arma::mat& W, 
                      float alpha, float beta);
double posterior_eta_cpp(const arma::vec& eta_n, const arma::ivec& v_n, const arma::mat& W, 
                       float alpha, float beta);

arma::vec grad_posterior_eta_cpp(const arma::vec& eta_n, const arma::ivec& v_n, 
                                 const arma::mat& W, const arma::rowvec& norms_W,
                                 float alpha, float beta);

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
  
#endif
