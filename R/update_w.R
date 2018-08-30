#' @title Maximize w.r.t W using H with multiplicative updates
#' @param W the (whole) current dictionary matrix.
#' @param H activation coefficients matrix (J,K,N)
#' @param V matrix to be factorized.
update_W_mc_H_depprecated <- function(W, V, H_samples){
  nsamples <- dim(H_samples)[1]
  numerator <- 0
  for(i in 1:nsamples){
    numerator <- numerator + (V/(W%*%H_samples[i,,]))%*%t(H_samples[i,,])
  }
  denominator <- diag(1/apply(H_samples, MARGIN=2, FUN=sum)) # sum over J,N 
  W * (numerator%*%denominator)
}

#' @title Maximize w.r.t W using H with multiplicative updates
#' @param W the current dictionary matrix.
#' @param H activation coefficients matrix (J,K,N)
#' @param V matrix to be factorized.
update_W_mc_H <- function(W, V, H_samples, updates=1){
  nsamples <- dim(H_samples)[1]
  
  for(u in 1:updates){
    numerator <- 0
    for(i in 1:nsamples){
      numerator <- numerator + (V/(W%*%H_samples[i,,]))%*%t(H_samples[i,,])
    }
    denominator <- diag(1/apply(H_samples, MARGIN=2, FUN=sum)) # sum over J,N 
    W <- W * (numerator%*%denominator)
  }
  W
}

#' @title Maximize w.r.t W using C,H
#' @param C_samples samples of C tensor (F,K,N)
#' @param H_samples samples of H matrix (K,N)
update_W_mc_CH <- function(C_samples_mean, H_samples_mean){
  diag_H_k <- diag(1/base::rowMeans(H_samples_mean)) #  mean over N
  base::rowMeans(C_samples_mean, dims=2) %*% diag_H_k # mean over N
}

#' @title Maximize w.r.t W given C
#' @param W the (whole) current dictionary matrix.
#' @param C_samples samples of the C tensor
#' @param alpha parameter alpha
#' @param beta parameter beta
update_W_mc_C <- function(C_samples_mean, alpha=1, beta=1){
  beta/alpha * base::rowMeans(C_samples_mean, dims = 2) # mean over N
}

#' @title Maximize w.r.t W using C,H
#' @param C_mean weighted means of C (F,K, niter)
#' @param H_samples samples of H matrix (K, niter)
update_W_mc_CH_SAEM <- function(C_weighted_sum, H_weighted_sum){
  C_weighted_sum %*% diag(1/H_weighted_sum)
}

#' @title Maximize w.r.t W using C,H
#' @param C_mean weighted means of C (F,K, niter)
#' @param H_samples samples of H matrix (K, niter)
update_W_mc_C_SAEM <- function(C_weighted_sum, alphabeta_weighted_sum){
  C_weighted_sum/alphabeta_weighted_sum
}