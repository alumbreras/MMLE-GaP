#' @title Hamiltonian sample of a column h_n
#' @param v_n n-th column of matrix V (F x 1)
#' @param W dictionary matrix (F x K)
#' @param h_n_current (K x 1) initial value h_n (n-th column of matrix H)
#' @param epsilon length of each Hamiltonian step
#' @param L number of Hamiltonian steps at each iteration ("path length")
#' @param iter number of samples (or "iterations")
sample_hamiltonian_h <-function(v_n, W, h_n_current, epsilon = 0.1, L=10, iter=100,
                                alpha, beta){
  if(length(v_n) != dim(W)[1])         stop("Non-conforming dimension F in v, W")
  if(length(h_n_current) != dim(W)[2]) stop("Non-conforming dimension K in h, W")
  
  norms_W <- apply(W, 2, sum) # pre-compute
  
  # # re-use the common posterior function to avoid stupid errors
  U <- function(h_n, v_n, W){
       #-posterior_h(h_n, v_n, W, alpha, beta, log = TRUE)
        -posterior_h_cpp(h_n, v_n, W, alpha, beta)
   }
  
  grad_U <- function(h_n, v_n, W, norms_W){
    K <-dim(W)[2]
   #res <-  -c((alpha-1)/h_n - 
   #       (beta + norms_W) +
   #        t(v_n %*% (W/((W%*%h_n)%*%rep(1,K)))))
    
    #rescpp <- grad_U_cpp(h_n, v_n, W, norms_W, alpha, beta)
    #cat("\nres: ", res)
    #cat("\nrescpp: ", rescpp)
    #cat("\nh_n: ", h_n)
    grad_U_cpp(h_n, v_n, W, norms_W, alpha, beta)
  }
  
  grad_U_R <- function(h_n, v_n, W, norms_W, alpha, beta){
    K <-dim(W)[2]
    res <-  -c((alpha-1)/h_n - 
           (beta + norms_W) +
            t(v_n %*% (W/((W%*%h_n)%*%rep(1,K)))))
    res
  }
  
  
  h_n_samples <- array(NA, dim=c(iter, length(h_n_current)))
  h_n_samples[1,] <- h_n_current # for transparency, init point at the beginning of the chain
  
  current_q <- h_n_current # we name it q since it is common practice
  for(i in 2:iter){
    q <- current_q
    p <- rnorm(length(q), 0,1 ) # independent standard normal variates
    current_p <- p
    
    
    # Make a half step for momentum at the beginning
    p <- p - epsilon * grad_U(q, v_n, W, norms_W) / 2

    # Alternate a full step for the position
    for(l in 1:L){
      
      # Make a full step for the position
      q <- q + epsilon * p
      q <- abs(q) # Bouncing. Recommended by Nico.
      
      # Make a full step for the momentum, except at the end of trajectory
      if(l != L) {
        p <- p - epsilon * grad_U(q, v_n, W, norms_W)
      }
    }
    # Make a half step for momentum at the end.
    p <- p - epsilon * grad_U(q, v_n, W, norms_W) / 2
    
    # Negate momentum at end of trajectory to make the proposal symmetric
    p <- -p
    
    # Evaluate potential and kinetic energies at start and end of trajectory
    current_U  <- U(current_q, v_n, W)
    current_K <- sum(current_p^2) / 2
    proposed_U <- U(q, v_n, W)
    proposed_K <- sum(p^2) / 2
    
    # Accept or reject the state at the end of trajectory, returning either
    # the position at the end of the trajectory or the initial position
    current_H <- current_U + current_K
    proposed_H <- proposed_U  + proposed_K
    
    if(runif(1) < exp(current_H -  proposed_H)){
      current_q <- q  # accept
    }
    h_n_samples[i,] <- current_q
  }
  return(h_n_samples)
}