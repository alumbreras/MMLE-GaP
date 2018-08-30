# Approximate computations of the marginal likelihood p(V|W)

#' @title Likelihood p(V | W) with Chib's method
#' @param V original matrix
#' @param W estimated dictionary matrix
#' @param H any H matrix
#' @param C samples of matrix of auxiliary variables c
#' @details Computes equation 22 of the paper Dikmen, Févotte (2012)
#' to compute p(H | V, W) we use the samples from the last iterations
#' since these are closer to the true posterior than the older ones where 
#' W was not optimal yet.
marginal_likelihood_Chib <- function(V, W, H, C, alpha, beta){
  F <- dim(V)[1]
  N <- dim(V)[2]
  K <- dim(W)[2]

  WH <- W%*%H
  likelihood <- 0
  
  # If K=1, H is coarced into a vector, but we want to keep its matrix form
  if(K == 1){
    H <- t(H)
  }
  
  # log p(V | W, H*)
  for(i in 1:F){
    for(j in 1:N){
      likelihood <- likelihood + dpois(V[i,j], WH[i,j], log = TRUE)
    }
  }
  # log p(H*)
  for(i in 1:K){
    for(j in 1:N){
      likelihood <- likelihood + dgamma(H[i,j], alpha, scale = beta, log = TRUE)
    }
  }
  
  # log p(H | V, W)
  likelihood_h_post <- 0
  nsamples <- dim(C)[4]  
  for(i in 1:nsamples){
	C_sum <- colSums(C[,,,i], dims=1)
	alpha_post <- alpha + C_sum # matrix of alpha_post (KN)
	beta_post  <- 1/(beta + colSums(W)) # vector of beta_post (K)
    for(k in 1:K){
      for(n in 1:N){
        likelihood_h_post <- likelihood_h_post + dgamma(H[k,n], 
                                                        alpha_post[k,n], 
                                                        scale = beta_post[k], 
                                                        log = FALSE)
      }
    }
  }
  likelihood <- likelihood - log(likelihood_h_post/nsamples)
  likelihood
}

#' @title Marginal likelihood with samples from the prior
#' @param V data matrix 
#' @param W dictionary matrix
#' @param alpha parameter for the prior
#' @param beta parameter for the prior
#' @param nsamples number of samples from the Importance distribution
marginal_likelihood_basic <- function(V, W, alpha=1, beta=1, nsamples=100){

  F <- dim(W)[1]  
  K <- dim(W)[2]
  N <- dim(V)[2]
  H_samples <- array(NA, dim = c(K, N, nsamples))
  
  # Draw samples from the prior
  for(n in 1:N){
    for(k in 1:K){
      H_samples[k,n,] <- rgamma(nsamples, shape = alpha, rate = beta)
    }
  }
  
  # Compute the approximation to p(V|W)
  logp <- 0
  for(n in 1:N){
    prob_n <- 0
    for(i in 1:nsamples){
       H <- H_samples[,,i]
       WH <- W %*% H
       prob_n <- prob_n + exp(sum(sapply(1:F, function(f) dpois(V[f,n], WH[f,n], log=TRUE))))
    }
    prob_n <- prob_n/nsamples
    logp <- logp + log(prob_n)
  }
  logp
}


library(MGLM)
#' @title Marginal likelihood with Importance Sampling
#' @param V data matrix 
#' @param W dictionary matrix
#' @param alpha parameter for the Negative Multinomial
#' @param beta parameter for the Negative Multinomial
#' @param nsamples number of samples from the Importance distribution
marginal_likelihood_is <- function(V, W, alpha=1, beta=1, nsamples=100){
  
  F <- dim(W)[1]  
  K <- dim(W)[2]
  N <- dim(V)[2]
  C_samples <- array(NA, dim = c(F, K, N, nsamples))
  
  # Draw samples from a importance distribution
  probs_multinomial <- t(apply(W, 1, function(x) x/sum(x)))
  for(n in 1:N){
    for(f in 1:F){
      C_samples[f,,n,] <- rmultinom(nsamples, size = V[f,n], prob = probs_multinomial[f,])
    }
  }
  
  # Compute the approximation to p(V|W)
  probs_nm <- apply(W, 2, function(x) x/(beta + sum(x)))
  logp <- 0
  for(n in 1:N){
    meanprob <- 0
    for(i in 1:nsamples){
      logp_nb <- sum(sapply(1:K, function(k) {
        dnegmn(C_samples[,k,n,i], probs_nm[,k], alpha)}
      )
      )
      logp_q <- sum(sapply(1:F, function(f) {
        dmultinom(C_samples[f,,n,i], prob = probs_multinomial[f,], log = TRUE)
      }
      )
      )
      meanprob <- meanprob + exp(logp_nb - logp_q - log(nsamples))
    }
    logp <- logp + log(meanprob) # sum logprob to the other N data points 
  }
  logp
}
