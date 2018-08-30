# This implements MCEM inference for the paper
# Dikmen, Févotte, 2012. Maximum Marginal Likelihood Estimation for 
# Nonnegative Dictionary Learning in the Gamma-Poisson Model. 
# In IEEE Transactions on Signal Processing, Vol. 60, NO. 10, October 2012.

#TODO: sample_z must return only post-burnin samples
#TODO: sample_metropolis must return only post-burnin samples

accepted <- 0
proposed  <- 0
accepted_ <- 0
proposed_  <- 0

#' @title NMF with MMLE
#' @param V matrix to be factorized
#' @param K number of latent factors (or dimensions)
#' @param W initial matrix W
#' @param maxiters maximum number of iterations
#' @param nsamples number of samples per iteration in the E-step
#' @param mc mcem or saem.
#' @details The number of iterations is set to maxiters.
#' In the future, there should be a convergence check in case
#' the KL converges before maxiters.
#' 
#' If saem is used, then the number of samples gradually increases until reaching
#' nsamples for the last iterations
#' @return W last estimator of W
#'         traces$H: a K x N x samples tensor
#'         traces$W: a K x samples - matrix vector with the evolution of row norms 
nmf_mmle_mcem <- function(V, initW, maxiters = 100, nsamples = 50, 
                          algo = "gibbs", mincrease="constant", alpha=1, beta=1,
                          mstep ="C", mult_updates = 1, saem=FALSE,
                          L=10, epsilon=0.5, burnin=0.5){
    
  if(algo == 'gibbs_par'){
    library(foreach)
    library(parallel)
    library(doParallel)
    cl <- makeCluster(4)
    registerDoParallel(cl)
  }
  
  F <- dim(V)[1]
  K <- dim(initW)[2]
  N <- dim(V)[2]
  m <- rep(nsamples, maxiters) # samples per iteration
  npostburnin <- floor(nsamples*burnin) # burnin samples
  
  
  if(dim(V)[1] != dim(initW)[1]) {
    stop("Dimensions between V and initW do not match")
  }

  if((nsamples<20) && (mincrease != "constant")){
    warning("Number of samples too small. Setting mincrease = 'constant'")
    mincrease = "constant"
  }
  
  if(nsamples<2) stop('Too few samples by iteration')
  

  # C = FKN nsamples because the C++ imposes it (arma::icube object has
  # the slices in the last dimension)
  C_samples_mean <- array(NA, dim=c(F, K, N))
  H_samples_mean <- array(1, dim=c(K, N)) 
  H_samples <- array(1, dim=c(npostburnin, K, N))
  W_traces <- array(NA, dim=c(maxiters, F, K))

  # Storage of means for SAEM updates
  C_means <- array(NA, dim=c(maxiters, F, K))
  H_means <- array(NA, dim=c(maxiters, K))
  alphabeta <- rep(NA, maxiters)
  
  times <- rep(NA, maxiters)

  W <- initW
  
  # Initialize the first matrix using row-wise multinomials
  # so that we have a valid matrix where sum_k C = V
  # TODO: slow in big data. 
  # for(f in 1:F){
  #   cat("\n f:", f)
  #   for(n in 1:N){
  #     #C_samples_mean[f,,n] <- t(rmultinom(1, size=V[f,n], prob=rep(1/K, K)))[1,]
  #     C_samples_mean[f,,n] <- V[f,n]/K
  #   }
  # }
  
  
  # If SAEM: start with few samples and gradually increase (E-step)
  if(saem){
    
    # Linear increse:
    # samples per iteration (linear increase)
    # grow such that total number of samples is equivalent
    # to the non-saem version
    totalsamples <- nsamples*maxiters
    minsamples <- 20
    maxsamples <- 2*totalsamples/maxiters
    m <- minsamples + (1:maxiters) * (maxsamples - minsamples)/maxiters
    m <- sapply(m, floor)
    
    # Constant
    m <- rep(nsamples, maxiters)
    
    # classic weights 
    if(mincrease == 'linear'){
      m <- minsamples+(1:maxiters)*(nsamples - minsamples)/maxiters
      m <- sapply(m, floor)
    } 
    
    k_saem <- floor(0.5*maxiters) # at 1.0 is EM
    weights <- c(rep(1,k_saem), 1/(((k_saem+1):maxiters)-k_saem))
    
    # exponential decay weights à la Dual-Averaging / NUTS
    # used here: https://www.tandfonline.com/doi/pdf/10.1198/106186006X157469
    alpha_saem <- 0.75
    weights <- (1/1:maxiters)^alpha_saem
  }
  
  accepted <- 0
  rejected <- 0
  niter <- 1
  idx_trace <- 1
  for (niter in 1:maxiters){
    
    cat("\n ***************")
    cat("\n * niter",niter, "/", maxiters)
    cat("\n ***************")

      
    # Gibbs --------------------------------------------------------------------
    if(algo == 'gibbs'){
      cat("\n E-step")
      # E-step
      for(n in 1:N){
          #if(!n %% 2){cat('\n iter:', niter, '. Gibbs in column ', n, "/", N)}
          if(!n %% 1000){cat('\ngibbs in column ', n, "/", N)}
          samples <- sample_gibbs_cpp(V[,n], W, H_samples_mean[,n],
                                      alpha=alpha, beta=beta,
                                      iter=m[niter], burnin=0.5)
          # H_samples[nsamples,,] is not compatible with SAEM where nsamples is variable
          # Thus, we only use this when strictly necessary
          if(mstep == "H" || mstep == "Hupdates") {
            H_samples[,,n] <- samples$h_n_samples
          }
          C_samples_mean[,,n] <- samples$C_n_samples_mean
          H_samples_mean[,n]  <- samples$h_n_samples_mean
          #chain <- mcmc(H_samples[,,n])
          #plot(chain)
      }
      
      # M-step
      cat("\n M-step")
      if(!saem){
        if(mstep == "H"){
          W <- update_W_mc_H(W, V, H_samples, updates = mult_updates)
        } else if(mstep == "CH") {
          W <- update_W_mc_CH(C_samples_mean, H_samples_mean)
        } else if(mstep == "C"){
          W <- update_W_mc_C(C_samples_mean, alpha, beta)
        } else {
          stop("Unknown mstep")
        }
      }
      if(saem){
        
        # Store E[C], E[H]... for current iteration
        C_means <- C_means*(1-weights[niter])
        C_means[niter,,] <- weights[niter] * rowMeans(C_samples_mean, dims=2)
        
        H_means <- H_means*(1-weights[niter])
        H_means[niter,] <- weights[niter] * rowMeans(H_samples_mean)
        
        alphabeta <- alphabeta*(1-weights[niter])
        alphabeta[niter] <- weights[niter] * alpha/beta
        
        C_weighted_sum <- colSums(C_means, na.rm = TRUE)
        H_weighted_sum <- colSums(H_means, na.rm = TRUE)
        alphabeta_weighted_sum <- sum(alphabeta, na.rm = TRUE)
        
        if(mstep == "CH") {
          W <- update_W_mc_CH_SAEM(C_weighted_sum, H_weighted_sum)
        } else if(mstep == "C"){
          W <- update_W_mc_C_SAEM(C_weighted_sum, alphabeta_weighted_sum)
        } else {
          stop("Unknown mstep for SAEM")
        } 
      }
    }
    
    # Gibbs parallel------------------------------------------------------------
    if(algo == 'gibbs_par'){
      # seems I need to load Rcpp functions in a different package, WTF
      # https://stackoverflow.com/questions/35090327/r-object-of-type-s4-is-not-subsettable-when-using-foreach-loop
      cat("\n E-step")
      # E-step
      parresults <- foreach( n = 1:N, .packages=c("gibbsMMLE", "Matrix")) %dopar% {
                        #sourceCpp("../src/sample_gibbs.cpp")
                        samples <- gibbsMMLE::sample_gibbs_cpp(V[,n], W, H_samples_mean[,n],
                                                     alpha=alpha, beta=beta,
                                                     iter=m[niter], burnin=0.5)
                        list(n=n, samples=samples)
      }
      for(n in 1:N){
        nn <- parresults[[n]]$n
        H_samples[,,nn]      <- parresults[[nn]]$samples$h_n_samples
        C_samples_mean[,,nn] <- parresults[[nn]]$samples$C_n_samples_mean
        H_samples_mean[,nn]  <- parresults[[nn]]$samples$h_n_samples_mean
      }
      
      # M-step
      cat("\n M-step")
      if(mstep == "H"){
        W <- update_W_mc_H(W, V, H_samples)
      } else if(mstep == "CH") {
        W <- update_W_mc_CH(C_samples_mean, H_samples_mean)
      } else if(mstep == "C"){
        W <- update_W_mc_C(C_samples_mean, alpha, beta)
      } else {
        stop("Unknown mstep")
      }
    }
    # Gibbs N --------------------------------------------------------------------
    if(algo == 'gibbs_n'){
      cat("\n E-step")
      # E-step
      samples <- sample_gibbs_cpp_N(V, W, H_samples_mean,
                                    alpha=alpha, beta=beta, 
                                    iter=m[niter], burnin=0.5)
      H_samples <- samples$H_samples
      C_samples_mean <- samples$C_samples_mean
      H_samples_mean  <- samples$H_samples_mean
      
      # M-step
      cat("\n M-step")
      if(mstep == "H"){
        W <- update_W_mc_H(W, V, H_samples)
      } else if(mstep == "CH") {
        W <- update_W_mc_CH(C_samples_mean, H_samples_mean)
      } else if(mstep == "C"){
        W <- update_W_mc_C(C_samples_mean, alpha, beta)
      } else {
        stop("Unknown mstep")
      }
    }
    
    # Gibbs Z--------------------------------------------------------------------
    if(algo == 'gibbs_z'){
      cat("\n E-step")
      # E-step
      for(n in 1:N){
        if(!n %% 2){cat('\n iter:', niter, '. Gibbs in column ', n, "/", N)}
        if(!n %% 1000){cat('\ngibbs in column ', n)}
        samples <- sample_gibbs_z_cpp(V[,n], W, C_samples_mean[,,n], 
                                    alpha=alpha, beta=beta, iter=m[niter],
                                    burnin=0.5)
        if(sum(is.na(samples))) stop("NA values in Gibbs sampler")
        C_samples_mean[,,n] <- samples
      }
      
      # M-step
      cat("\n M-step")
      W <- update_W_mc_C(C_samples_mean, alpha, beta)
    }
    
    
    # MH ----------------------------------------------------------------------
    if(algo == 'metropolis'){
      
      # E-step

        # sample H | W, V
        for(n in 1:N){
          samples <- sample_metropolis_h_reparam_cpp(V[,n], W, 
                                                  H_samples_mean[,n], 
                                                  alpha=alpha, 
                                                  beta=beta,
                                                  step=0.6, iter=m[niter]
          )
          H_samples[,,n] <- samples
          H_samples_mean[,n] <- base::colMeans(H_samples[,,n]) # over J samples
          
        }

      # M-step
      W <- update_W_mc_H(W, V, H_samples)
    }
    
    # HMC ---------------------------------------------------------
    if(algo =="hmc"){
      # E-step
      cat("\n E-step")
      for(n in 1:N){
        samples <- sample_hmc_cpp(V[,n], W, H_samples_mean[,n], 
                                  alpha=alpha, 
                                  beta=beta,
                                  L = L,
                                  epsilon= epsilon, iter=m[niter])
        H_samples[,,n] <- samples
        H_samples_mean[,n] <- base::colMeans(H_samples[,,n]) # over J samples
      }


      # M-step
      cat("\n M-step")
      W <- update_W_mc_H(W, V, H_samples)
    }
  
    # Hamiltonian Stan ---------------------------------------------------------
    if(algo == 'hamiltonian-stan'){
      stan_model <- rstan::stan_model("MMLE_H.stan")
      
      
      accepted <- 0
      total <- 0
      
      # E-step
      nsample <- 1      
      for(n in 1:N){
        # sample_h_n | W, V
        fit <- rstan::sampling(stan_model, 
                               data=list(K=K, F=F, W=W, v=V[,n]), 
                               iter=m[niter], chains=1, verbose=FALSE, 
                               algorithm = "NUTS", 
                               init=list(list(h=H_samples_mean[,n])))
        
        samples <- as.data.frame(rstan::extract(fit, permuted = FALSE, inc_warmup=FALSE)[,1,1:length(H_samples[1,,n])])
        H_samples[,,n] <- as.matrix(samples)
        H_samples_mean[,n] <- base::colMeans(H_samples[,,n]) # over J samples
        
      }

      # M-step
      cat("\n M-step")
      W <- update_W_mc_H(W, V, H_samples)
    }  
    
    W_traces[niter,,] <- W
    times[niter] <- Sys.time()
    

    } # end iters
    
  if(algo == 'gibbs_par'){
    stopCluster(cl)
  }
  
  
  times <- times - min(times)
  list(V=V, W=W, W_traces = W_traces, times = times)
}
