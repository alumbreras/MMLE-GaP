#' @title sample h_{kn} given all the other h in the same column 
#' with Metroplis Hastings
#' @param v_n n-th column of matrix V (F x 1)
#' @param W dictionary matrix (F x K)
#' @param h_n_current (K x 1) initial value h_n (n-th column of matrix H)
#' @param alpha parameter alpha
#' @param beta parameter beta
#' @param step step size for the proposal distribution
#' @param iter number of samples (or "iterations")
sample_metropolis_gaussian_h <- function(v_n, W, h_n_current, alpha=1, beta=1, step=1, iter=100){
  
  h_n_samples <- array(1, dim=c(iter, length(h_n_current)))
  h_n_samples[1,] <- h_n_current # for transparency, init point at the beginning of the chain
  
  for(i in 1:iter){
    h_n <- h_n_current
    for(k in 1:K){
      h_n[k] <- h_n[k] + rnorm(1, sd=step)
      a <- posterior_h(h_n, v_n, W, alpha, beta, log=TRUE) - 
           posterior_h(h_n_current, v_n, W, alpha, beta, log=TRUE)
      if(runif(1) < exp(a)){
        h_n_current <- h_n
      }
    }
    h_n_samples[i,] <- h_n_current
  }
  return(h_n_samples)
}


#' @title sample h_{n} given all the other h in the same column 
#' with Metroplis Hastings. Use a Gamma proposal.
#' @param v_n n-th column of matrix V (F x 1)
#' @param W dictionary matrix (F x K)
#' @param h_n_current (K x 1) initial value h_n (n-th column of matrix H)
#' @param alpha parameter alpha
#' @param beta parameter beta
#' @param step step size for the proposal distribution
#' @param iter number of samples (or "iterations")
sample_metropolis_gamma_h <- function(v_n, W, h_n_current, alpha=1, beta=1, step=1, iter=100){
  
  h_n_samples <- array(1, dim=c(iter, length(h_n_current)))
  
  transition_ratio <- 0
  for(i in 1:iter){
    h_n <- h_n_current
    for(k in 1:K){
      h_n[k] <- rgamma(1, shape= step, rate = step/h_n[k]) #E = shape/rate
      cat("\n", h_n_current[k], step,  h_n[k], step/h_n[k])
      transition_ratio <- transition_ratio + 
        dgamma(h_n_current[k], shape = step, rate = step/h_n[k], log=TRUE) -
        dgamma(h_n[k],         shape = step, rate = step/h_n_current[k], log=TRUE)
    }
  
    a <- posterior_h(h_n, v_n, W, alpha, beta, log=TRUE) - 
         posterior_h(h_n_current, v_n, W, alpha, beta, log=TRUE) +
         transition_ratio
  
    if(runif(1) < exp(a)){
      h_n_current <- h_n
    }
    h_n_samples[i,] <- h_n_current
  }
  return(h_n_samples)
}


#' @title sample h_{n} given all the other h in the same column 
#' with Metroplis Hastings. Use a log-transform reparamatrization
#' @param v_n n-th column of matrix V (F x 1)
#' @param W dictionary matrix (F x K)
#' @param h_n_current (K x 1) initial value h_n (n-th column of matrix H)
#' @param alpha parameter alpha
#' @param beta parameter beta
#' @param step step size for the proposal distribution
#' @param iter number of samples (or "iterations")
sample_metropolis_h_reparam <- function(v_n, W, h_n_current, alpha=1, beta=1, step=1, iter=100){
  
  posterior_eta <- function(eta_n, v_n, W, alpha, beta, log=TRUE){
    posterior_h(exp(eta_n), v_n, W, alpha, beta, log=TRUE) + 
      log(det(diag(exp(eta_n))))
  }
  
  eta_n_samples <- array(1, dim=c(iter, length(h_n_current)))
  eta_n_current <- log(h_n_current)
  K <- length(h_n_current)
  
  for(i in 1:iter){
    eta_n <- eta_n_current
    eta_n <- eta_n + rnorm(K, sd=step)
    #a <- posterior_h(exp(eta_n), v_n, W, alpha, beta, log=TRUE) + log(det(diag(exp(eta_n)))) - 
    #     posterior_h(exp(eta_n_current), v_n, W, alpha, beta, log=TRUE) - log(det(diag(exp(eta_n_current))))

    a <- posterior_eta(eta_n, v_n, W, alpha, beta, log=TRUE) -
         posterior_eta(eta_n_current, v_n, W, alpha, beta, log=TRUE)
    
    if(runif(1) < exp(a)){
      eta_n_current <- eta_n
    }
    eta_n_samples[i,] <- eta_n_current
  }
  return(exp(eta_n_samples))
}


# This method is very inefficient. Sampling H is better.
sample_metropolis_c <- function(nsamples=1, C_n, W, alpha=1, beta=1){
  
  F <- dim(C_n)[1]
  K <- dim(C_n)[2]
  
  C_samples <- array(NA, dim=c(nsamples,F, K))
  
  for(nsample in 1:nsamples){
    
    # Jump from last sample, or from the initial point
    # if first sample
    if(nsample == 1){
      C_last <- C_n 
    } else{
      C_last <- C_samples[nsample-1 ,,] 
    }
    C_new <- C_last
    
    # Propose new matrix C, modifying feature by feature
    # to guarantee the constraints (sum over K dimensions)
    for(f in 1:F){
      
      # if the vector is empty, the is no possible transition
      if(sum(C_last[f,]) == 0){
        next
      }
      
      # choose among the non-zero positions
      # If only the n-th position is non-zero,
      # sample(which(C_new[f,]>0), 1) choses from 1:n instead of choosing n
      # these two lines are the standard alternative
      allowed_pos <- which(C_new[f,]>0)
      k_minus <- allowed_pos[sample(length(allowed_pos), 1)]
      C_new[f,k_minus] <-  C_new[f,k_minus] - 1
      
      
      # choose among any position
      k_plus <- sample(K, 1)
      C_new[f,k_plus] <- C_new[f,k_plus] + 1
    }
    
    # Accept or reject
    prob_jump <- 0
    prob_jump_rev <- 0
    for(f in 1:F){
      prob_jump <- prob_jump - sum(C_n[f,]>0) - length(C_n[f,])
      prob_jump_rev <- prob_jump_rev - sum(C_new[f,]>0) - length(C_new[f,])
    }
    
    a <- prob_jump_rev - 
      prob_jump +
      likelihood_C_W(C_new,  W, alpha=1, beta=1, log=TRUE) -
      likelihood_C_W(C_last,  W, alpha=1, beta=1, log=TRUE)
    
    proposed <- proposed + 1
    if(runif(1) < exp(a)){
      C_sample <- C_new
      accepted <- accepted + 1
    } else{
      C_sample <- C_last
    }
    C_samples[nsample,,] <- C_sample
  }
  
  return(C_samples)
}