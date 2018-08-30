#' @title posterior of a column of H
#' @param h a column vector from H, of length K
#' @description proportionality constants are ignored. Thus, the posterior is
#' unnormalized. 
#' @details For the sake of simplicity, this is done by multiplying the gamma prior
#' and the poisson likelihood
posterior_h <- function(h, v, W, alpha=1, beta=1, log = FALSE){
  
  if(any(h<0)){
    logp <- -Inf
  }
  else{
    logp <- sum(dpois( v, lambda = W %*% h, log=TRUE)) + 
            sum(dgamma(h, shape = alpha, rate = beta, log=TRUE))
  }
  
  if(log == TRUE){
    return(logp)
  } else{
    return(exp(logp))
  }
}

posterior_h_1 <- function(h, v, W, alpha=1, beta=1, log = FALSE){
  
  if(any(h<0)){
    return(-log(0))
  }
  logp <- c(log(prod(h)^(alpha-1)) -
          h %*% (beta + apply(W, 2, sum)) +
          v %*% log(W%*%h))
  
  if(log == TRUE){
    return(logp)
  } else{
    return(exp(logp))
  }
}

posterior_h_2 <- function(h, v, W, alpha=1, beta=1, log = FALSE){
  
  if(any(h<0)){
    logp <- -Inf
  }
  else{
    F <- dim(W)[1]
    K <- dim(W)[2]
    logp <- 0
    for(f in 1:F){
      logp <- logp + dpois(v[f],  lambda = sum(W[f,] * h), log=TRUE)
    }
    for(k in 1:K){
      logp <- logp + dgamma(h[k],  shape = alpha, rate = beta, log=TRUE)
    }
  }
  
  if(log == TRUE){
    return(logp)
  } else{
    return(exp(logp))
  }
}


