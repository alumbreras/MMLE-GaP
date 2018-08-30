# Compute the value p(v | W) using the 
# expression of the close form
###############################################
library(partitions)
library(foreach)
library(doParallel)
library(abind)

#' @title Negative Multinomial distribution with weight parameters
#' @param x vector of counts
#' @param w vector of weights in every dimension
#' @param alpha number of failure experiments
#' @param beta weight of the failure dimension.
dnegmul_prints <- function(x, w=1:length(x), alpha=1, beta=1, log=FALSE){
  prob <- 0
  denominator <- beta + sum(w)
  cat('\n denominator:', denominator)
  cat('\n C_norm:', sum(x))
  
  prob <- prob + lgamma(alpha + sum(x)) - lgamma(alpha) - sum(lgamma(x+1))
  cat('\n Gamma:', lgamma(alpha + sum(x)) - lgamma(alpha) - sum(lgamma(x+1)), '\n')
  
  prob <- prob + log((beta/denominator)^alpha)
  cat('\n Beta:', log((beta/denominator)^alpha), '\n')
  
  prob <- prob + sum(log((w/denominator)^x))
  cat('\n W:',  sum(log((w/denominator)^x)), '\n')
  
  if(log == FALSE){
    prob <- exp(prob)
  }
  prob
}

dnegmul <- function(x, w=1:length(x), alpha=1, beta=1, log=FALSE){
  prob <- 0
  denominator <- beta + sum(w)
  prob <- prob + lgamma(alpha + sum(x)) - lgamma(alpha) - sum(lgamma(x+1))
  prob <- prob + log((beta/denominator)^alpha)
  prob <- prob + sum(log((w/denominator)^x))
  if(log == FALSE){
    prob <- exp(prob)
  }
  prob
}

#' @title Negative Multinomial distribution with probability parameters
dnegmul_prob <- function(x, p=1:length(x), alpha=1, log=FALSE){
  prob <- 0
  prob <- prob + lgamma(alpha + sum(x)) - lgamma(alpha) - sum(lgamma(x+1))
  prob <- prob + log((1-sum(p))^alpha)
  prob <- prob + sum(log(p^x))
  if(log == FALSE){
    prob <- exp(prob)
  }
  prob
}


#' @title P(c | W) when alpha is infinite (F=1)
likelihood_c_W_inf <- function(c, lambda){
  exp(-sum(lambda)) * prod(lambda^c / prod(factorial(c)))
}

#' @title P(v | W) when alpha is infinite (F=1)
likelihood_v_W_inf <- function(v, lambda){
  partitions <- compositions(n=v , m=2)
  prob <- 0
  for (i in 1:ncol(partitions)){
    prob <- prob  +  prod(lambda^partitions[,i]) / prod(factorial(partitions[,i]))
  }
  prob <- exp(-sum(lambda)) * prob
  prob
}

#' @title P(v | W) when alpha is infinite (F=1)
likelihood_V_W_inf <- function(V, lambda){
  prob <- 0
  for (n in 1:length(V)){
    prob <- prob + log(likelihood_v_W_inf(V[n], lambda))
  }
  exp(prob)
}

#' @title Product of independent Negative Multinomial Distributions
#' @param C matrix of counts
#' @param W matrix of weights in every dimension
#' @param alpha number of failure experiments
#' @param beta 1/weight of the failure dimension.
#' @details Column k of C is associated to column k of W.
likelihood_C_W <- function(C, W, alpha=1, beta=1, log=FALSE){
  K <- ncol(W)
  prob <- 0 
  for (k in 1:K){
    prob <- prob + dnegmul(C[,k], W[,k], alpha, beta, log=TRUE)
  }

  if(log == FALSE){
    prob <- exp(prob)
  }
  prob
}

#' @title likelihood p(v | W)
#' @param v a column of V where every component is a feature.
#' @description computes p(v|W) using the non-simplified expression
#' @return log probability
#' @details Computes the likelihood of a column of V. Since each column
#' is independent, the likelihood of the whole matrix V is the product of the
#' likelihoods of each column in V.
likelihood_v_W <- function(v, W, alpha, beta, log=TRUE, vectorizable=TRUE){
  F <- length(v)
  K <- ncol(W)
  
  # Precompute all partitions for every feature
  partitions <- list()
  for(f in 1:F){
    partitions[[f]] <- compositions(n=v[f],m=K)#[:,1:2]
  }
  npartitions <- sapply(partitions, ncol)
  N <- prod(npartitions)
  indices <- lapply(npartitions, function(x) 1:x)
  grid <- as.matrix(do.call(expand.grid, indices)) # breaks if too big

  # Pick one partition from each feature and compute its probability
  prob <- 0
  for(n in 1:nrow(grid)){
      if(n %% 10000 == 0) cat('\n partition:', n)
    
      selected <- c(grid[n,])
      C <- t(sapply(1:F, function(f) partitions[[f]][,selected[f]]))
      prob <- prob + likelihood_C_W(C, W, alpha, beta, log = FALSE)

  }
  if(log==TRUE){
    prob <- log(prob)  
  }
  res <- list(prob = prob,
              iterations = N)
  if(vectorizable){
    res <- prob
  }
  res
}

#' @title likelihood p(v | W)
#' @param v a column of V where every component is a feature.
#' @description computes p(v|W) using the simplified expression (if alpha=1)
#' and F = 1
#' @return log probability
#' @details Computes the likelihood of a column of V. Since each column
#' is independent, the likelihood of the whole matrix V is the product of the
#' likelihoods of each column in V.
likelihood_v_W_alpha1 <- function(v, W, beta, log=TRUE, vectorizable=TRUE){
  F <- length(v)
  K <- ncol(W)
  denominators <- beta + colSums(W)
  
  prob <- 0
  for(f in 1:F){
    partitions <- compositions(n=v[f],m=K)
  
    prob_partition <- 0
    for(n in 1:ncol(partitions)){
      if(n %% 10000 == 0) cat('\n partition:', n)
      c <- partitions[,n]
    
      prob_k <- 0
      for (k in 1:K){
        prob_k <- prob_k + (1/F) * (log(beta) - log(denominators[k]))
        #prob_k <- prob_k + c[k] * (log(W[f,k]) - log(denominators[k]))
        prob_k <- prob_k + log((W[f,k])^c[k]) - c[k] * log(denominators[k])
        
      }
      prob_partition <- prob_partition + exp(prob_k)
    }
    prob <- prob + log(prob_partition)
  }
  prob <- exp(prob)

  if(log==TRUE){
    prob <- log(prob)  
  }
  res <- list(prob = prob,
              iterations = N)
  if(vectorizable){
    res <- prob
  }
  res
}


#' @title likelihood p(V | W)
#' @param V a full matrix with one observation per column
#' @description computes p(V|W)
#' @return log probability
#' @details Computes the likelihood of a full V. Since each column
#' is independent, the likelihood of the whole matrix V is the product of the
#' likelihoods of each column in V.
likelihood_V_W <- function(V, W, alpha, beta, log=TRUE, vectorizable = FALSE){
  N <- ncol(V)
  F <- nrow(W)
  prob <- 0
  
  if((alpha == 1) && (F == 1)){
    # Faster computation
    for (n in 1:N){
      prob <- prob + likelihood_v_W_alpha1(V[,n], W, beta, log=TRUE, vectorizable=TRUE)
    }
  }
  else{
    for (n in 1:N){
      prob <- prob + likelihood_v_W(V[,n], W, alpha, beta, log=TRUE, vectorizable=TRUE)
    }
  }
  
  
  if(log == FALSE){
    prob <- exp(prob)
  }
  prob
    
}

#' @title manual computation for F=1, K=2, N=2, and e.g. v=[1,6]
#' @details this is the most simple example that shows two modes
#' The goal of this function is to "debug" the model and understand 
#' at what point the bimodality emerges.
#' @param C: c vectors to be used. C is N x K dimensions
likelihood_C_W_manual <- function(C, w, alpha=10){
  
  prob <- 1
  
  prob <- prob * 
          gamma(alpha + C[1,1]) * gamma(alpha + C[1,2]) * 
          gamma(alpha + C[2,1]) * gamma(alpha + C[2,2])
  prob <- prob / (gamma(alpha) * prod(factorial(C)))
  
  prob <- prob * (1 + w[1])^(-2*alpha - colSums(C)[1])
  prob <- prob * (1 + w[2])^(-2*alpha - colSums(C)[2])
  prob <- prob * w[1]^colSums(C)[1] * w[2]^colSums(C)[2]
  prob
}

#' @title C partitions generator
#' @param v A feature vector ( a column of V)
#' @details Compute all possible matrices C that explain a feature vector
#' @return An array that represents a cube. Each slice in the cube is a matrix C.
compute_partitions <- function(v, K = 2){
  #v <- t(t(c(1,2)))
  F <- nrow(v)
  partitions <- list()
  for(f in 1:F){
    # all vectors c explaining this feature
    partitions[[f]] <- compositions(n=v[f],m=K)
  }
  npartitions <- sapply(partitions, ncol)
  indices <- lapply(npartitions, function(x) 1:x)
  grid <- as.matrix(do.call(expand.grid, indices))
  N <- prod(npartitions)
  #N <- nrow(grid)
  C_cube <- array(dim=c(F,K,N))
  for(n in 1:nrow(grid)){
    selected <- c(grid[n,])
    C_cube[,,n] <- t(sapply(1:F, function(f) partitions[[f]][,selected[f]]))
  }
  
  C_cube
}

#' title: Compute C partitions for the equation
#' where every term has all N
compute_partitions_N <- function(V, K = 2){

  F <- nrow(V)
  C_cubes <- list()
  cubes <- list()
  N <- ncol(V)
  for(n in 1:N){
    cubes[[n]] <- compute_partitions(V[,n, drop=FALSE], K = K)
  }
  npartitions <- sapply(cubes, function(x) dim(x)[3])
  indices <- lapply(npartitions, function(x) 1:x)
  grid <- as.matrix(do.call(expand.grid, indices))

  for(nworld in 1:nrow(grid)){
    selected <- c(grid[nworld,])
    res <- lapply(1:N, function(n) cubes[[n]][,,selected[n]])
    cube <- abind(res, rev.along=0)
    if(F==1){
      cube <- t(cube) # N x K
    }
    else{
      cube <- aperm(cube, c(1,3,2))
    }
    C_cubes[[nworld]] <- cube
  }
  C_cubes
}

compute_partitions_N_strengths <- function(V, K=2){
  #V <- matrix(c(1,2), nrow=1)
  partitions <- compute_partitions_N(V, K)
  partitions <- abind(partitions, rev.along=0)
  c_weights <- t(apply(partitions, 3, function(x) colSums(x)))
  c_weights[order(c_weights[,1]),]
}

likelihood_V_W_manual <- function(V, w, alpha=10){
  partitions <- compute_partitions_N(V, K=2)
  
  # Pick one partition from each feature and compute its probability
  prob <- 0
  for(n in 1:length(partitions)){
    C <- partitions[[n]]
    prob <- prob + likelihood_C_W_manual(C, w, alpha=10)
  }
  prob
}

#' p(C | W) using sum over N and F=1
likelihood_C_W_2 <- function(C, p, alpha=1){
  
  prob <- 1
  K <- dim(C)[1]
  N <- dim(C)[2]
  
  for(n in 1:N){
    for(k in 1:K){
      prob <- prob * 
        gamma(alpha + C[k,n]) * gamma(alpha + C[k,n]) / 
       ( gamma(alpha) * prod(factorial(C[k,n])) )
    }
  }
  
  for(k in 1:K){
    prob <- prob * (1 - p[k])^(N*alpha) * p[k]^(rowSums(C)[k])
  }
  prob
}

#' Likelihood p(V | W) using sum over N and F=1
likelihood_V_W_2 <- function(V, p, alpha=10){
  partitions <- compute_partitions_N(V, K=2)
  
  # Pick one partition from each feature and compute its probability
  prob <- 0
  for(n in 1:length(partitions)){
    C <- t(partitions[[n]]) # transpose to make it K x N
    prob <- prob + likelihood_C_W_2(C, p, alpha=1)
  }
  prob
}




TEST <- FALSE
if(TEST){
  # Data
  ################################################
  # Synthetic 
  usertopics <- rep(1:5, each=2)
  movietopics <- rep(1:5, each=10)
  W <- matrix(0, 10, 5)
  H <- matrix(0, 5, 50)
  for(i in 1:10){
    W[i,usertopics[i]] <- 5  
  }
  for(i in 1:50){
    H[movietopics[i],i] <- 5
  }
  V <- W%*%H
  V <- ceiling(V/2)
  
  # Louis
  V <- as.matrix(read.csv(file = 'dataset_1.csv', header=FALSE))
  colnames(V) <- NULL
  
  ################################################
  alpha <- 1
  beta <- 1
  K <- 5
  W <- matrix(1, nrow=nrow(V), ncol=K)
  
  like <- 0
  iters <- 0
  #for (n in 1:ncol(V)){
  for (n in 1:100){
    res <- likelihood_v_W(V[,n], W, alpha, beta, log=TRUE, vectorizable = FALSE)
    like <- like + res$prob
    cat('\n n:', n, ": ", res$prob, "///",  like, "...")
    iters <- iters + res$iterations
  }
  
  #profvis(likelihood_V_W(c(20,2), W, alpha, beta))
  cat("\nProbability: ", like)
  cat("\nIterations:", iters)

  likelihood_V_W(V[,1:100], W, alpha=1, beta, log=TRUE, vectorizable = FALSE)
  likelihood_V_W(V, W, alpha=1, beta, log=TRUE, vectorizable = FALSE)
  }