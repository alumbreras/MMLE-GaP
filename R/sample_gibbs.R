#' @title Conditional sampler of c_fn
#' @param v_fn value of V[f,n] the drawn vector will sum up to this value
#' @param w_f f-th row of the W matrix
#' @param h_n h-th column of the H matrix
#' @description Samples a multinomial vector c_fn
#' @details Deals with cases where v_fn = 0, therefore c_fn = 0,..0
sample_gibbs_c <- function(v_fn, w_f, h_n){
  K <- length(h_n)
  probs <- w_f * h_n
  probs_denominator <- sum(probs)
  if(sum(probs_denominator)==0){return(rep(0, K))}
  c(rmultinom(1, v_fn, probs))
}

#' @title Conditional sampler of h_kn
#' @param w_k  k-th column of matrix W
#' @param c_kn vector in position (k,n) of tensor C [F x K x N]
sample_gibbs_h <- function(w_k, c_kn, alpha=1, beta=1){
  alpha_post <- alpha + sum(c_kn)
  beta_post <- 1/(1/beta + sum(w_k))
  rgamma(1, alpha_post, scale=beta_post)
}

#' @title Gibbs sampler of h_n using auxiliary variables C
#' @param v_n n-th column of matrix V (F x 1)
#' @param W dictionary matrix (F x K)
#' @param h_n (K x 1) initial value h_n (n-th column of matrix H)
#' @param alpha parameter alpha
#' @param beta parameter beta
#' @param iter number of samples (or "iterations")
sample_gibbs <- function(v_n, W, h_n, alpha=1, beta=1, iter=100){
  
  if(dim(W)[2] != length(h_n)) stop ("K lengths do not match")
  
  K <- length(h_n)
  F <- length(v_n)
  C_n_samples <- array(NA, dim=c(F, K, iter))
  h_n_samples <- array(1, dim=c(iter, K))
  h_n_samples[1,] <- h_n  # for transparency, init point at the beginning of the chain
  C_n <- array(NA, dim=c(F, K))
  h_n <- rep(1, K)
  
  for(i in 2:iter){
    for(f in 1:F){
      C_n[f,] <- sample_gibbs_c(v_n[f], W[f,], h_n)
    }
    for(k in 1:K){
      h_n[k] <- sample_gibbs_h(W[,k], C_n[,k], alpha, beta)
    }
    C_n_samples[,,i] <- C_n
    h_n_samples[i,] <- h_n
  }
  return(list(C_n_samples = C_n_samples, h_n_samples = h_n_samples))
}