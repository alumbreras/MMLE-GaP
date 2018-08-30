# Given a non-negative matrix V, find a decomposition  V = WH using the 
# multiplicative updates given in :
# Lee, D. D., & Seung, H. S. (2000). Algorithms for Non-negative Matrix Factorization. 
# In Advances in Neural Information Processing Systems (NIPS’2000) (Vol. 13, pp. 556–562).

#' @title Update W_fk
#' @param f row of W
#' @param k column of W
#' @param V matrix being factorized
#' @param W dictionary matrix
#' @param H activation coefficients matrix 
#' @details Updates the given position using Lee and Seung multiplicative updates
update_w_Lee <- function(f, k, V, W, H){
  numerator <- sum(H[k,]*V[f,]/(W%*%H)[f,], na.rm = TRUE)
  denominator <- sum(H[k,])
  if(denominator == 0) stop("H matrix with zero-sum column")
  W[f,k] <- W[f,k] * numerator/denominator
  W[f,k]
}

#' @title Update H_kn
#' @param k row of H
#' @param n column of H
#' @param V matrix being factorized
#' @param W dictionary matrix
#' @param H activation coefficients matrix 
#' @details Updates the given position using Lee and Seung multiplicative updates
update_h_Lee <- function(k, n, V, W, H){
  numerator <- sum(W[,k]*V[,n]/(W%*%H)[,n], na.rm = TRUE)
  denominator <- sum(W[,k])
  if(denominator == 0) stop("W matrix with zero-sum column")
  H[k,n] <- H[k,n] * numerator/denominator
  H[k,n]
}

#' @title Compute the KL divergence 
#' @param A first matrix
#' @param B second matrix
#' @details This function is useful for debugging
KL_divergence_Lee <- function(A, B){
  I <- dim(V)[1]
  J <- dim(V)[2]
  total <- 0
  for(i in 1:I){
    for(j in 1:J){
      if (A[i,j] == 0) next
      if ((B[i,j]) == 0 && (A[i,j] > 0)) return(Inf)
      total <- total + A[i,j]*log(A[i,j]/B[i,j]) - A[i,j] + B[i,j]
    }
  }
  total
}

#' @title NMF with Lee and Saung multiplicative updates
#' @param V matrix to factorize
#' @param K number of latent factors (or dimensions)
#' @param W initial W matrix
#' @param H initial H matrix
#' @param maxiters maximum number of iterations
#' @details The number of iterations is set to maxiters.
#' In the future, there should be a convergence check in case
#' the KL converges before maxiters
nmf_Lee <- function(V, K, W, H, maxiters = 100){
  KLlog <- rep(NA, maxiters)
  for (i in 1:maxiters){
    # update H
    for(k in 1:K){
      for(n in 1:N){
        H[k,n] <- update_h_Lee(k, n, V, W, H)
      }
    }
    # update W
    for(f in 1:F){
      for(k in 1:K){
        W[f,k] <- update_w_Lee(f, k, V, W, H)
      }
    }
    KLlog[[i]] <- KL_divergence_Lee(V, W%*%H)
  }
  list(W=W, H=H, KLlog = KLlog)
}


