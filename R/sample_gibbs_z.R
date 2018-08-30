#' @title Gibbs sampler of h_n using auxiliary variables C
#' @param v_n n-th column of matrix V (F x 1)
#' @param W dictionary matrix (F x K)
#' @param C_n (F x K) initial value C_n
#' @param alpha parameter alpha
#' @param beta parameter beta
#' @param iter number of samples (or "iterations")
sample_gibbs_z <- function(v_n, W, C_n, alpha=1, beta=1, iter=100){
  
  F <- dim(W)[1]
  K <- dim(W)[2]
  
  C_n_samples <- array(NA, dim=c(F, K, iter))
  probs <- rep(NA, K)
  denominators = beta + colSums(W);
  
  C_n_samples[,,1] = C_n;
  L = colSums(C_n);
  
  # Resample the topic assigment k of each count
  for(i in 2:iter){
    
    # Fix the original matrix so that we can sample each occurrence
    C_n_initial <- C_n
    
    for(f in 1:F){
      # v_n[f] total counts of feature f
      # for each occurent of f in the document and topic k, re-sample its topic
      for(k in 1:K){
        
        # skip empty counts
        if(C_n_initial[f,k] == 0) { 
          next
        }
      
        # Use a C_initial snapshot so that each occurrence is re-sampled only once
        for(l in 1:C_n_initial[f,k]){
            L[k] = L[k] - 1             # remove occurrence
            C_n[f, k] = C_n[f, k] - 1   # remove occurrence
            
            # Sample new k for the occurrence
            probs = (alpha + L + 1) * (W[f,]/denominators);
            probs = probs/sum(probs);
            k_chosen = sample(K, 1, replace=FALSE, prob=probs);
            
            C_n[f, k_chosen] = C_n[f, k_chosen] + 1   # add occurrence
            L[k_chosen] = L[k_chosen] + 1             # add occurrence
          }
        }
      }
      C_n_samples[,,i] = C_n;
    }
  return(C_n_samples)
}