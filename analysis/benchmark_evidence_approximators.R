# Compares how the differents estimators of  p(V|W) convergence to the true
# value when we increase the number of used samples
library(dplyr)
library(ggplot2)

# Generate synthetic data from a GaP model ------------------------------------
alpha <- 1
beta <- 1
F <- 3
K <- 2
N <- 100

W <- t(rdirichlet(K, alpha = rep(1,F)))*F
H <- matrix(rgamma(K*N, shape=alpha, rate=beta), ncol=N) # K x N matrix
V <- matrix(NA, nrow=F, ncol=N)
for(f in 1:F){
  for(n in 1:N){
    V[f,n] <- rpois(1, lambda = W[f,] %*% H[,n])
  }
}

truelike <-  likelihood_V_W(V, W, alpha, beta)

nsamples_list <- seq(10, 2000, by=100)
df <- data.frame()
for(nsamples in nsamples_list){
  cat("\nComputing with nsamples = ", nsamples)
  basic <- marginal_likelihood_basic(V, W, alpha, beta, nsamples)
  IS <- marginal_likelihood_is(V, W, alpha, beta, nsamples)
  df <- rbind(df, list(nsamples=nsamples, basic=basic-truelike, IS=IS-truelike))
}

df <- df %>%  gather(method, error, basic, IS)
ggplot(df, aes(x=nsamples, y=abs(error), group=method, color=method)) + 
  geom_point() + geom_line() +
  theme_bw()