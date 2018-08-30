
Q1 <- function(theta, E_C, E_lambdas){
  -theta*sum(E_lambdas) + log(theta) * sum(E_C)
}

Q2 <- function(theta, E_C){
  -(alpha + sum(E_C))*log(beta + theta) + sum(E_C)*log(theta)
}


gibbs <- function(x, theta, lambdas, alpha=1, beta=1){
  nsamples <- 100
  K <- nrow(lambdas)
  lambdas_samples <- array(NA, dim = c(nsamples, K, N))
  C_samples <- array(NA, dim = c(nsamples, K, N))
  C <- array(NA, dim = c(K, N))
  
  for(n in 1:N){
    #cat("\n N:", n)
    for(i in 1:nsamples){
      
      # sample c_n
      probs <- theta*lambdas[,n]
      C[,n] <- rmultinom(1, x[n], probs/sum(probs))
      
      # sample lambda_kn
      for(k in 1:K){
        alpha_k <- alpha + C[k,n]
        beta_k <- beta + theta
        lambdas[k,n] <- rgamma(1, alpha_k, rate = beta_k)
      }
      
      # save samples
      C_samples[i,,n] <- C[,n]
      lambdas_samples[i,,n] <- lambdas[,n]
    }
  }
  
  # remove burnin
  nburnin <- floor(nsamples/2)
  C_samples <- C_samples[-c(1:nburnin),,] 
  lambdas_samples <- lambdas_samples[-c(1:nburnin),,]
  
  return(list(E_C = colMeans(C_samples, dims=1), 
              E_lambdas = colMeans(lambdas_samples, dims=1)))
}

EM_CH <- function(x, K, theta, iter, alpha=1, beta=1){
  
  V <- t(as.matrix(x))
  thetas <- seq(0.1,50, by=0.1)
  y <- sapply(thetas, function(theta) {likelihood_V_W_cpp(V, as.matrix(theta), alpha, beta)})
  df <- data.frame(theta=thetas, value = y, method="true")
  plot(thetas, y)
  N <- length(x)
  E_lambdas <- array(1, dim = c(K,N))
  
  
  for(i in 1:iter){
    cat("\nem iter:", i)
    samples <- gibbs(x, theta, E_lambdas, alpha, beta)
    E_C <- samples$E_C
    E_lambdas <- samples$E_lambdas
    theta = sum(E_C)/sum(E_lambdas)
    cat("\ntheta:", theta)
    
    lb <- Q1(thetas, E_C, E_lambdas)
    df <- rbind(df, data.frame(theta=thetas, value = lb, method=as.character(i)))
    p <- ggplot(df, aes(x = theta, y=value, group=method, color=method)) + geom_line() + theme_bw() +
          xlab("theta") + ylab("Q-CH") + ggtitle(i)
    print(p)
    plot(thetas, y-lb)
  }
}

EM_C <- function(x, K, theta, iter, alpha=1, beta=1){
  
  V <- t(as.matrix(x))
  thetas <- seq(0.1,50, by=0.1)
  y <- sapply(thetas, function(theta) {likelihood_V_W_cpp(V, as.matrix(theta), alpha, beta)})
  df <- data.frame(theta=thetas, value = y, method="true")
  
  N <- length(x)
  E_lambdas <- array(1, dim = c(K,N))
  for(i in 1:iter){
    cat("\nem iter:", i)
    cat("\nE_lambdas dim", dim(E_lambdas))
    samples <- gibbs(x, theta, E_lambdas, alpha, beta)
    E_C <- samples$E_C
    E_lambdas <- samples$E_lambdas
    theta = (beta/alpha) * sum(E_C)/N
    theta = (beta/alpha) * mean(x) # if w_1=...=w_k (does not depend on the Gibbs sampler)
    
    cat("\ntheta:", theta)
    
    lb <- Q2(thetas, E_C)
    df <- rbind(df, data.frame(theta=thetas, value = lb, method=as.character(i)))
    p <- ggplot(df, aes(x = theta, y=value, group=method, color=method)) + geom_line() + theme_bw() +
      xlab("theta") + ylab("Q-C") + ggtitle(i)
    print(p)
  }
}

EM_H <- function(x, K, theta, iter, alpha=1, beta=1){
  
  V <- t(as.matrix(x))
  thetas <- seq(0.1,50, by=0.1)
  y <- sapply(thetas, function(theta) {likelihood_V_W_cpp(V, as.matrix(theta), alpha, beta)})
  df <- data.frame(theta=thetas, value = y, method="true")
  
  N <- length(x)
  E_lambdas <- array(1, dim = c(K,N))
  for(i in 1:iter){
    cat("\nem iter:", i)
    cat("\nE_lambdas dim", dim(E_lambdas))
    samples <- gibbs(x, theta, E_lambdas, alpha, beta)
    E_C <- samples$E_C
    E_lambdas <- samples$E_lambdas
    theta = mean(x) / sum(E_lambdas)
    cat("\ntheta:", theta)
    
    lb <- Q2(thetas, E_C)
    df <- rbind(df, data.frame(theta=thetas, value = lb, method=as.character(i)))
    p <- ggplot(df, aes(x = theta, y=value, group=method, color=method)) + geom_line() + theme_bw() +
      xlab("theta") + ylab("Q-C") + ggtitle(i)
    print(p)
  }
}

# Generate data from Gamma-Poisson mixture
N <- 1000
theta_star <- 65
alpha <- 10
lambdas <- rgamma(N, alpha, 1)
x <- rpois(N, lambda = theta_star * lambdas)
hist(x)
cat("mean:", mean(x))
cat("var: ", var(x))


likelihood_V_W_cpp(V, W, alpha, beta)



# Estimate the parameter theta
EM_CH(x, K=5, theta=1, iter=35, alpha, beta)