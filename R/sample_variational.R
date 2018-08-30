# Variational inference of p(H, C | V, W)

# This implements Variational Bayes EM inference for the paper
# Dikmen, Févotte, 2012. Maximum Marginal Likelihood Estimation for 
# Nonnegative Dictionary Learning in the Gamma-Poisson Model. 
# In IEEE Transactions on Signal Processing, Vol. 60, NO. 10, October 2012.
#
# _var v_stands por variational parameters
#alpha <- 1
#beta <- 1



lower_bound <- function(V, W, E_c, E_ln_h, E_h, q_c_probs, q_h_alpha, q_h_beta,
                        alpha, beta){
  F <- dim(V)[1]
  N <- dim(V)[2]
  K <- dim(W)[2]
  
  pC <- 0
  pH <- 0
  entropy_C <- 0
  entropy_H <- 0 
  
  for(n in 1:N){
    # terms from p(C)
    for (f in 1:F){
      # na.rm to avoid 0*-Inf and set it directly to zero
      pC <- pC - 
            W[f,] %*% E_h[,n] + 
            sum(E_c[f,n,] * ( log(W[f,]) + E_ln_h[,n]), na.rm=TRUE)
    }
    
    # terms from p(H)
    for(k in 1:K){
      pH <- pH + 
            (alpha-1) * E_ln_h[k,n] - 
            E_h[k,n] * beta + 
            alpha * log(beta) - 
            lgamma(alpha)
    }
    
    # terms from Entropy(C)
    for(f in 1:F){
      entropy_C <- entropy_C - 
                 lfactorial(V[f,n]) -  
                 sum(E_c[f,n,] * (log(q_c_probs[f,n,])), na.rm=TRUE)
    }
    
    # terms from Entropy(H)
    for(k in 1:K){
      var_alpha <- q_h_alpha[k,n]
      var_beta  <- q_h_beta[k,n]
      entropy_H <- entropy_H - 
              lgamma(var_alpha) + 
              (var_alpha - 1) * digamma(var_alpha) + 
              log(var_beta) - var_alpha
    }
  }
  lower <- pC + pH - entropy_C - entropy_H
  lower
}


# Uses shape/rate  parametrization (E = alpha/beta)
optim_variational <- function(V, W, maxiters = 100, alpha = 1, beta = 1){
  
  update_q_c <- function(f, n, W, E_ln_h){
    K <- dim(W)[2]
    var_probs <- W[f,] * exp(E_ln_h[,n])
    var_probs <- var_probs/sum(var_probs)
    list(var_probs = var_probs)
  }
  
  update_q_h <- function(k, n, W, E_c, alpha, beta){
    var_alpha <- alpha + sum(E_c[,n,k])
    var_beta <- beta + sum(W[,k])
    list(var_alpha = var_alpha, var_beta = var_beta)
  }
  
  
  lower_bound_log <- rep(NA, maxiters)
  V <- as.matrix(V, nrow = dim(W)[1])
  F <- dim(V)[1]
  N <- dim(V)[2]
  K <- dim(W)[2]
  
  # Structures to store the estimated parameters
  q_c_probs <- array(NA, dim=c(F, N, K))
  q_h_alpha <- array(NA, dim=c(K, N))
  q_h_beta <- array(NA, dim=c(K, N))
  
  E_c <- array(1/K, dim=c(F, N, K))
  E_h <- array(runif(K*N), dim=c(K, N)) # random init
  E_ln_h <- array(runif(K*N), dim=c(K, N)) # random init
  W <- array(runif(F*K), dim=c(F, K)) # random init
  
  for (iter in 1:maxiters){
    cat("\niteration", iter)
    # Update variational posterior q_c 
    for(n in 1:N){
      for(f in 1:F){
        q_c_probs[f,n,] <- update_q_c(f, n, W, E_ln_h)$var_probs
        E_c[f,n,] <- q_c_probs[f,n,] * V[f,n]
      }
    }
    
    # update variational posterior q_h
    for(n in 1:N){
      for(k in 1:K){
        params <- update_q_h(k, n, W, E_c, alpha, beta)
        q_h_alpha[k,n] <- params$var_alpha
        q_h_beta[k,n]  <- params$var_beta
        E_h[k,n] <- q_h_alpha[k,n] / q_h_beta[k,n]
        E_ln_h[k,n] <- digamma(q_h_alpha[k,n]) - log(q_h_beta[k,n])
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
  list(lower_bound = lower_bound_log, q_h_alpha=q_h_alpha, q_h_beta=q_h_beta)
  
}


