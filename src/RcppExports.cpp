// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// posterior_h_cpp
double posterior_h_cpp(const arma::vec& h_n, const arma::ivec& v_n, const arma::mat& W, float alpha, float beta);
RcppExport SEXP _rMMLE_posterior_h_cpp(SEXP h_nSEXP, SEXP v_nSEXP, SEXP WSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type h_n(h_nSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_h_cpp(h_n, v_n, W, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// posterior_eta_cpp
double posterior_eta_cpp(const arma::vec& eta_n, const arma::ivec& v_n, const arma::mat& W, float alpha, float beta);
RcppExport SEXP _rMMLE_posterior_eta_cpp(SEXP eta_nSEXP, SEXP v_nSEXP, SEXP WSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta_n(eta_nSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_eta_cpp(eta_n, v_n, W, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// grad_posterior_eta_cpp_R
SEXP grad_posterior_eta_cpp_R(const arma::vec& eta_n, const arma::ivec& v_n, const arma::mat& W, const arma::rowvec& norms_W, float alpha, float beta);
RcppExport SEXP _rMMLE_grad_posterior_eta_cpp_R(SEXP eta_nSEXP, SEXP v_nSEXP, SEXP WSEXP, SEXP norms_WSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type eta_n(eta_nSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type norms_W(norms_WSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_posterior_eta_cpp_R(eta_n, v_n, W, norms_W, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP _rMMLE_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// likelihood_V_W_cpp
double likelihood_V_W_cpp(const arma::imat& V, const arma::mat& W, float alpha, float beta);
RcppExport SEXP _rMMLE_likelihood_V_W_cpp(SEXP VSEXP, SEXP WSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::imat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_V_W_cpp(V, W, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// sample_gibbs_cpp
SEXP sample_gibbs_cpp(const arma::vec& v_n, const arma::mat& W, arma::vec h_n, double alpha, double beta, int iter, double burnin);
RcppExport SEXP _rMMLE_sample_gibbs_cpp(SEXP v_nSEXP, SEXP WSEXP, SEXP h_nSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP iterSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h_n(h_nSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_gibbs_cpp(v_n, W, h_n, alpha, beta, iter, burnin));
    return rcpp_result_gen;
END_RCPP
}
// sample_gibbs_cpp_N
SEXP sample_gibbs_cpp_N(const arma::sp_mat& V, const arma::mat& W, arma::mat H_init, double alpha, double beta, int iter, double burnin);
RcppExport SEXP _rMMLE_sample_gibbs_cpp_N(SEXP VSEXP, SEXP WSEXP, SEXP H_initSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP iterSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H_init(H_initSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_gibbs_cpp_N(V, W, H_init, alpha, beta, iter, burnin));
    return rcpp_result_gen;
END_RCPP
}
// sample_gibbs_z_cpp
SEXP sample_gibbs_z_cpp(const arma::vec& v_n, const arma::mat& W, arma::imat C, double alpha, double beta, int iter, double burnin);
RcppExport SEXP _rMMLE_sample_gibbs_z_cpp(SEXP v_nSEXP, SEXP WSEXP, SEXP CSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP iterSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::imat >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type burnin(burninSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_gibbs_z_cpp(v_n, W, C, alpha, beta, iter, burnin));
    return rcpp_result_gen;
END_RCPP
}
// sample_hmc_cpp
SEXP sample_hmc_cpp(const arma::ivec& v_n, const arma::mat& W, arma::vec h_n_current, float alpha, float beta, int L, float epsilon, int iter);
RcppExport SEXP _rMMLE_sample_hmc_cpp(SEXP v_nSEXP, SEXP WSEXP, SEXP h_n_currentSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP LSEXP, SEXP epsilonSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec& >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h_n_current(h_n_currentSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< float >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_hmc_cpp(v_n, W, h_n_current, alpha, beta, L, epsilon, iter));
    return rcpp_result_gen;
END_RCPP
}
// sample_mala_cpp
SEXP sample_mala_cpp(const arma::ivec v_n, const arma::mat W, const arma::vec h_n_current, double alpha, double beta, double delta, int iter);
RcppExport SEXP _rMMLE_sample_mala_cpp(SEXP v_nSEXP, SEXP WSEXP, SEXP h_n_currentSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP deltaSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::ivec >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type h_n_current(h_n_currentSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_mala_cpp(v_n, W, h_n_current, alpha, beta, delta, iter));
    return rcpp_result_gen;
END_RCPP
}
// sample_metropolis_h_reparam_cpp
SEXP sample_metropolis_h_reparam_cpp(arma::ivec v_n, const arma::mat& W, arma::vec h_n_current, float alpha, float beta, float step, int iter);
RcppExport SEXP _rMMLE_sample_metropolis_h_reparam_cpp(SEXP v_nSEXP, SEXP WSEXP, SEXP h_n_currentSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP stepSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h_n_current(h_n_currentSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type step(stepSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_metropolis_h_reparam_cpp(v_n, W, h_n_current, alpha, beta, step, iter));
    return rcpp_result_gen;
END_RCPP
}
// posterior_h_cpp_
double posterior_h_cpp_(arma::vec& h_n, const arma::ivec& v_n, const arma::mat& W, float alpha, float beta);
RcppExport SEXP _rMMLE_posterior_h_cpp_(SEXP h_nSEXP, SEXP v_nSEXP, SEXP WSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type h_n(h_nSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< float >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(posterior_h_cpp_(h_n, v_n, W, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}
// sample_nuts_cpp
SEXP sample_nuts_cpp(arma::ivec v_n, const arma::mat& W, arma::vec h_n_current, double alpha, double beta, float epsilon, int iter);
RcppExport SEXP _rMMLE_sample_nuts_cpp(SEXP v_nSEXP, SEXP WSEXP, SEXP h_n_currentSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP epsilonSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type h_n_current(h_n_currentSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_nuts_cpp(v_n, W, h_n_current, alpha, beta, epsilon, iter));
    return rcpp_result_gen;
END_RCPP
}
// variational_ch_cpp
SEXP variational_ch_cpp(arma::ivec v_n, arma::mat W, arma::vec alpha_var, arma::vec beta_var, double alpha, double beta, int maxiters);
RcppExport SEXP _rMMLE_variational_ch_cpp(SEXP v_nSEXP, SEXP WSEXP, SEXP alpha_varSEXP, SEXP beta_varSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP maxitersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::ivec >::type v_n(v_nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alpha_var(alpha_varSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_var(beta_varSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiters(maxitersSEXP);
    rcpp_result_gen = Rcpp::wrap(variational_ch_cpp(v_n, W, alpha_var, beta_var, alpha, beta, maxiters));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rMMLE_posterior_h_cpp", (DL_FUNC) &_rMMLE_posterior_h_cpp, 5},
    {"_rMMLE_posterior_eta_cpp", (DL_FUNC) &_rMMLE_posterior_eta_cpp, 5},
    {"_rMMLE_grad_posterior_eta_cpp_R", (DL_FUNC) &_rMMLE_grad_posterior_eta_cpp_R, 6},
    {"_rMMLE_mvrnormArma", (DL_FUNC) &_rMMLE_mvrnormArma, 3},
    {"_rMMLE_likelihood_V_W_cpp", (DL_FUNC) &_rMMLE_likelihood_V_W_cpp, 4},
    {"_rMMLE_sample_gibbs_cpp", (DL_FUNC) &_rMMLE_sample_gibbs_cpp, 7},
    {"_rMMLE_sample_gibbs_cpp_N", (DL_FUNC) &_rMMLE_sample_gibbs_cpp_N, 7},
    {"_rMMLE_sample_gibbs_z_cpp", (DL_FUNC) &_rMMLE_sample_gibbs_z_cpp, 7},
    {"_rMMLE_sample_hmc_cpp", (DL_FUNC) &_rMMLE_sample_hmc_cpp, 8},
    {"_rMMLE_sample_mala_cpp", (DL_FUNC) &_rMMLE_sample_mala_cpp, 7},
    {"_rMMLE_sample_metropolis_h_reparam_cpp", (DL_FUNC) &_rMMLE_sample_metropolis_h_reparam_cpp, 7},
    {"_rMMLE_posterior_h_cpp_", (DL_FUNC) &_rMMLE_posterior_h_cpp_, 5},
    {"_rMMLE_sample_nuts_cpp", (DL_FUNC) &_rMMLE_sample_nuts_cpp, 7},
    {"_rMMLE_variational_ch_cpp", (DL_FUNC) &_rMMLE_variational_ch_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_rMMLE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
