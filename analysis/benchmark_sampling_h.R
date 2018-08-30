# Benchmark multiple Monte Carlo samplers for a column pf H 
# using all the dictionary and the corresponding column of V.
library(dplyr)
library(tidyr)
library(fields)
library(rgl)
library(gtools)
library(rstan)
library(coda)

nsamples <- 1000

# ******************************************************************************

# Data 2 -----------------------------------------------------------------------
alpha <- 1
beta <- 1
K <- 2
N <- 10
W <- matrix(c(1,0,
              0,1), 
            ncol=K) # F x K matrix

# Data 1 -----------------------------------------------------------------------
alpha <- 1
beta <- 1
K <- 7
N <- 10
W <- matrix(c( 1,1,0,0,0,0,0,0,0,0,0,0,0,0,
               0,0,1,1,0,0,0,0,0,0,0,0,0,0,
               0,0,0,0,1,1,0,0,0,0,0,0,0,0,
               0,0,0,0,0,0,1,1,0,0,0,0,0,0,
               0,0,0,0,0,0,0,0,1,1,0,0,0,0,
               0,0,0,0,0,0,0,0,0,0,1,1,0,0,
               0,0,0,0,0,0,0,0,0,0,0,0,1,1), 
                   ncol=K) # F x K matrix
W <- W*100

# Data 3 -----------------------------------------------------------------------
alpha <- 1
beta <- 1
K <- 10
N <- 10
F <- 100
alphadir <- 1 # high: uniform vectors # low (<1): sparse vectors
                # inference easier over W with ortogonal rows
W <- rdirichlet(F, alpha = rep(alphadir,K))
W <- W*100
# ******************************************************************************

# Generate H from a Gamma and V = Poisson(WH)
h <- rgamma(K, shape=alpha, rate=1)
v <- rpois(F, lambda = W%*%h)


# initial sample for h
h_init <- rgamma(K, shape=alpha, rate=1) # K-length initialization of h column sampler

# ******************************************************************************
all.traces <- list()
times <- list()
# Metropolis -------------------------------------------------------------------
start <- Sys.time()
current_q <- h_init
traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=1, beta=1, step=0.01, iter=nsamples)
all.traces$metropolis <- traces
end <- Sys.time()
times$metropolis <- end-start

# Hamiltonian -----------------------------------------------------------------
# start <- Sys.time()
# L <- 10
# epsilon <- 0.008
# current_q <- h_init
# traces <- array(NA, dim=c(nsamples, length(current_q)))
# for(nsample in 1:nsamples){
#   current_q <- sample_hamiltonian_h(v, W, current_q, epsilon, L)
#   cat("\n", current_q)
#   traces[nsample,] <- current_q
# }
# all.traces$hamiltonian <- traces
# end <- Sys.time()
# times$hamiltonian <- end-start

# Hamiltonian NUTS STAN---------------------------------------------------------
start <- Sys.time()
model <- stan_model("MMLE_H.stan")
fit <- sampling(model, data=list(K=dim(W)[2], F=dim(W)[1], W=W, v=v), iter=nsamples, chains=1, 
                verbose=TRUE, algorithm = "NUTS")
traces <- as.data.frame(rstan::extract(fit)$h)
all.traces$nuts <- traces
end <- Sys.time()
times$nuts <- end-start


# Gibbs cpp
start <- Sys.time()
traces <- sample_gibbs_cpp_scalable(v, W, h_init, alpha=1, beta=1, iter=nsamples)$h_n_samples
all.traces$gibbscpp <- traces
end <- Sys.time()
times$gibbscpp <- end-start



# Plots ------------------------------------------------------------------------
for(i in 1:length(all.traces)){
  method <- names(all.traces)[i]
  cat("\nmethod:", method)
  
  # Plot traces
  traces <- all.traces[i]
  df.ground_truth <- data.frame(k=1:length(h), h=h)
  df <- as.data.frame(traces)
  df$sample <- 1:nrow(df)
  df <- gather(df, variable, value, -sample)
  df <- df[complete.cases(df),]
  p <- ggplot(df, aes(x=sample, y=value, color=variable)) + geom_line() + 
    geom_hline(data=df.ground_truth, aes(yintercept=h)) + theme_bw()
  print(p)
}

# Plot Autocorrelations
if(FALSE){
all.ess <- list()
for(i in 1:length(all.traces)){
  method <- names(all.traces)[i]
  cat("\nmethod:", method)
  
  # Plot autocorrelations
  traces <- all.traces[[i]]
  mcmc.traces <- mcmc(traces)
  autocorr.plot(mcmc.traces, lag.max = 100, auto.layout = TRUE)
  ess <- effectiveSize(mcmc.traces)
  print(ess)
  all.ess[[method]] <- mean(ess)
}

# Compare traces: autocorrelation, variance, CPU time, samples to converge
res <- list()
for(method in names(all.traces)){
  likelihoods <- apply(all.traces[[method]], 1, function(x) posterior_H(x, v, W, alpha, beta, log=TRUE))
  res[[method]] <- data.frame(method = method, sample=1:length(likelihoods), likelihood = likelihoods )
}
df.likelihoods <- do.call("rbind", res)
p <- ggplot(df.likelihoods, aes(x=sample, y=likelihood, group=method, color=method)) + 
  geom_line() + theme_bw()
print(p)


# Compare Squared Error (quite complementary to plot likelihoods)
res <- list()
for(method in names(all.traces)){
  errors <- apply(all.traces[[method]], 1, function(x) sqrt(t(x-h)%*%(x-h)))
  res[[method]] <- data.frame(method = method, sample=1:length(likelihoods), error = errors )
}
df.errors <- do.call("rbind", res)
p <- ggplot(df.errors, aes(x=sample, y=error, group=method, color=method)) + 
  geom_line() + theme_bw()
print(p)
}