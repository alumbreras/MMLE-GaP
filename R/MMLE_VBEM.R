# This implements Variational Bayes EM inference for the paper
# Dikmen, Févotte, 2012. Maximum Marginal Likelihood Estimation for 
# Nonnegative Dictionary Learning in the Gamma-Poisson Model. 
# In IEEE Transactions on Signal Processing, Vol. 60, NO. 10, October 2012.
#
# _var v_stands por variational parameters
alpha <- 1
beta <- 1

update_q_c <- function(f, n, W, E_ln_h){
  K <- dim(W)[2]
  var_probs <- W[f,] * exp(E_ln_h[,n])
  var_probs <- var_probs/sum(var_probs)
  list(var_probs = var_probs)
}

update_q_h <- function(k, n, W, E_c, alpha, beta){
  var_alpha <- alpha + sum(E_c[,n,k])
  var_beta <- 1/(sum(W[,k]) + 1/beta)
  list(var_alpha = var_alpha, var_beta = var_beta)
}

update_w_vb <- function(f, k, E_c, E_h){
  w <- sum(E_c[f,,k])/sum(E_h[k,])
}


lower_bound <- function(V, W, E_c, E_ln_h, E_h, q_c_probs, q_h_alpha, q_h_beta,
                        alpha, beta){
  F <- dim(V)[1]
  N <- dim(V)[2]
  K <- dim(W)[2]
  
  pC <- 0
  pH <- 0
  entropy <- 0 
  
  for(n in 1:N){
    # terms from p(C)
    for (f in 1:F){
      pC <- pC - W[f,] %*% E_h[,n]
      # na.rm to avoid 0*-Inf and set it directly to zero
      pC <- pC + sum(E_c[f,n,] * ( log(W[f,]) + E_ln_h[,n]), na.rm=TRUE)
    }
    
    # terms from p(H)
    for(k in 1:K){
      pH <- pH + (alpha-1) * E_ln_h[k,n]
      pH <- pH - E_h[k,n] / beta 
      pH <- pH - alpha * log(beta) - lgamma(alpha)
    }

    # terms from Entropy(C)
    # note there is a c! that cancelled out
    for(f in 1:F){
      var_probs  <- q_c_probs[f,n,]
      entropy <- entropy - lfactorial(V[f,n]) -  sum(E_c[f,n,] * (log(var_probs)), na.rm=TRUE)
    }
    
    # terms from Entropy(H)
    for(k in 1:K){
      var_alpha <- q_h_alpha[k,n]
      var_beta  <- q_h_beta[k,n]
      entropy <- entropy + lgamma(var_alpha)
      entropy <- entropy - (var_alpha - 1) * digamma(var_alpha)
      entropy <- entropy + log(var_beta)
      entropy <- entropy + var_alpha
    }
  }
  lower <- pC + pH + entropy
  lower
}


nmf_mmle_vbem <- function(V, K=5, maxiters = 100, alpha = 1, beta = 1){

  lower_bound_log <- rep(NA, maxiters)
  
  # Structures to store the estimated parameters
  F <- dim(V)[1]
  N <- dim(V)[2]
  
  q_c_probs <- array(NA, dim=c(F, N, K))
  q_h_alpha <- array(NA, dim=c(K, N))
  q_h_beta <- array(NA, dim=c(K, N))
  
  E_c <- array(1/K, dim=c(F, N, K))
  E_h <- array(runif(K*N), dim=c(K, N)) # random init
  E_ln_h <- array(runif(K*N), dim=c(K, N)) # random init
  W <- array(runif(F*K), dim=c(F, K)) # random init
  
  for (iter in 1:maxiters){
    # Update variational posterior q_c 
    for(f in 1:F){
      for(n in 1:N){
        q_c_probs[f,n,] <- update_q_c(f, n, W, E_ln_h)$var_probs
        E_c[f,n,] <- q_c_probs[f,n,] * V[f,n]
      }
    }
    
    # update variational posterior q_h
    for(k in 1:K){
      for(n in 1:N){
        params <- update_q_h(k, n, W, E_c, alpha, beta)
        q_h_alpha[k,n] <- params$var_alpha
        q_h_beta[k,n]  <- params$var_beta
        E_h[k,n] <- q_h_alpha[k,n] * q_h_beta[k,n]
        E_ln_h[k,n] <- digamma(q_h_alpha[k,n]) + log(q_h_beta[k,n])
      }
    }
    
    # update parameters W
    for(f in 1:F){
      for(k in 1:K){
        W[f,k] <- update_w_vb(f, k, E_c, E_h)
      }
    }
    
    # Monitor lower bound
    lower_bound_log[iter] <- lower_bound(V, W, E_c, E_ln_h, E_h, q_c_probs, 
                                         q_h_alpha, q_h_beta, alpha, beta)
  }
  if(sum(diff(lower_bound_log)<0)) {
    plot(lower_bound_log, main = "Shameful Lower Bound")
    stop("Sorry, that is embarrassing. Lower bound decreases at some point!")
  }
  list(lower_bound = lower_bound_log, W=W)
  
  
}