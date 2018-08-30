# Maximum Marginal Likelihood for Dictionary Learning in Gamma-Poisson models (MMLE-GaP)

Code related to the paper:

**Closed-form Marginal Likelihood in Gamma-Poisson Matrix Factorization**
Filstroff L., Lumbreras A., Févotte C. 
International Conference on Machine Learning (2018) 

R/Rcpp implementation of MMLE-GaP and related algorithms.

The main MCMC-EM algorithm is in R/MMLE_MCEM.R. From this file, calls are made to methods in other files that compute the Monte Carlo E-Step and the M-step. The basic MC E-Step is a Gibbs sampler, but we have also played with other alternatives. There are R and Rcpp implementations (src/ folder) that are much faster. 
