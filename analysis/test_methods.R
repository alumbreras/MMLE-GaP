# **************************************************************
# Script that launches all the implemented methods with toy data
# Use this to debug your implementations and to be sure everything
# is ok
# ***************************************************************

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(profvis)

set.seed(1)
maxiters <- 100
nsamples <- 200
alpha 	 <- 1
beta 	 <- 1
results <- list()

# Generate synthetic data from a GaP model ####
F <- 10
K <- 5
N <- 100
W <- matrix(c(10,10,0,0,0,0,0,0,0,0,
              0,0,10,10,0,0,0,0,0,0,
              0,0,0,0,10,10,0,0,0,0,
              0,0,0,0,0,0,10,10,0,0,
              0,0,0,0,0,0,0,0,10,10), 
            nrow=F,  ncol=K) # F x K matrix
F <- nrow(W)
K <- ncol(W)
H <- matrix(rgamma(K*N, shape=alpha, rate=1), ncol=N) # K x N matrix
V <- matrix(NA, nrow=F, ncol=N)
for(f in 1:F){
  for(n in 1:N){
    V[f,n] <- rpois(1, lambda = W[f,] %*% H[,n])
  }
}
plot_multiple_W(V, W, H, cols=3)

initW <- matrix(1, nrow = F, ncol = K+5) # K initially overestimated

# VARIATIONAL INFERENCE ------------------------------------------------
res_vb <- nmf_mmle_vbem(V, 10, maxiters = 500, alpha=alpha, beta=1)
plot_multiple_W(initW, res_vb$W, W)
par(mfrow=c(1,1))
plot(sort(colSums(res_vb$W), decreasing = TRUE)) # plot column norms

# METROPOLIS -----------------------------------------------------------
start <- Sys.time()
V <- V[,colSums(V)!=0] # drop empty columns
res_mh <- nmf_mmle_mcem(V, initW, maxiters, nsamples, 
                        algo='metropolis', mincrease = "constant", alpha=alpha)
end <- Sys.time()
cat("\nEllapsed time: ", end - start)

plot_multiple_W(initW, res_mh$W, W)
plot(sort(colSums(res_mh$W), decreasing = TRUE)) # plot column norms

results[['MH']] <- as.data.frame(apply(res_mh$W_traces, c(1,3), sum)) %>%
  mutate(mstep = "MH", iteration=1:maxiters, time=res_mh$times) %>%
  gather(k, value, -mstep, -iteration, -time)

# GIBBS C ----------------------------------------------------------------
start <- Sys.time()
res_gibbs <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep = "C")
end <- Sys.time()
cat("\nEllapsed time: ", end - start)

plot_multiple_W(initW, res_gibbs$W, W)
plot(sort(colSums(res_gibbs$W), decreasing = TRUE)) # plot column norms

results[["Gibbs C"]] <- as.data.frame(apply(res_gibbs$W_traces, c(1,3), sum)) %>%
  mutate(mstep = "Gibbs C", iteration=1:maxiters, time=res_gibbs$times) %>%
  gather(k, value, -mstep, -iteration, -time)

# GIBBS CH ----------------------------------------------------------------
start <- Sys.time()
res_gibbs <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep = "CH")
end <- Sys.time()
cat("\nEllapsed time: ", end - start)

plot_multiple_W(initW, res_gibbs$W, W)
plot(sort(colSums(res_gibbs$W), decreasing = TRUE)) # plot column norms

results[["Gibbs CH"]] <- as.data.frame(apply(res_gibbs$W_traces, c(1,3), sum)) %>%
  mutate(mstep = "Gibbs CH", iteration=1:maxiters, time=res_gibbs$times) %>%
  gather(k, value, -mstep, -iteration, -time)

# GIBBS H ----------------------------------------------------------------
start <- Sys.time()
res_gibbs <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep = "H")
end <- Sys.time()
cat("\nEllapsed time: ", end - start)

plot_multiple_W(initW, res_gibbs$W, W)
plot(sort(colSums(res_gibbs$W), decreasing = TRUE)) # plot column norms

results[["Gibbs H"]] <- as.data.frame(apply(res_gibbs$W_traces, c(1,3), sum)) %>%
  mutate(mstep = "Gibbs H", iteration=1:maxiters, time=res_gibbs$times) %>%
  gather(k, value, -mstep, -iteration, -time)


# HAMILTONIAN ----------------------------------------------------------
start <- Sys.time()

#profvis({
res_hmc <- nmf_mmle_mcem(V, initW, maxiters, nsamples, 
                         algo='hmc', 
                         mincrease = "constant", 
                         alpha=alpha, mstep ="H", L=10, epsilon=0.01)
#})
end <- Sys.time()
cat("\nEllapsed time: ", end - start)

plot_multiple_W(initW, res_hmc$W, W)
plot(sort(colSums(res_hmc$W), decreasing = TRUE)) # plot column norms

results[["HMC"]] <- as.data.frame(apply(res_hmc$W_traces, c(1,3), sum)) %>%
  mutate(mstep = "HMC", iteration=1:maxiters, time=res_hmc$times) %>%
  gather(k, value, -mstep, -iteration, -time)

# NUTS STAN---------------------------------------------------------------------
start <- Sys.time()

#profvis({
res_nuts <- nmf_mmle_mcem(V, initW, maxiters, nsamples, 
                         algo='hamiltonian-stan', 
                         mincrease = "constant", 
                         alpha=alpha, mstep ="H", L=10, epsilon=0.01)
#})
end <- Sys.time()
cat("\nEllapsed time: ", end - start)

plot_multiple_W(initW, res_hmc$W, W)
plot(sort(colSums(res_nuts$W), decreasing = TRUE)) # plot column norms

results[["HMC-STAN"]] <- as.data.frame(apply(res_nuts$W_traces, c(1,3), sum)) %>%
  mutate(mstep = "NUTS-STAN", iteration=1:maxiters, time=res_nuts$times) %>%
  gather(k, value, -mstep, -iteration, -time)


################################################################################
# Plot traces of column norms (for all methods)
################################################################################
df <- do.call(rbind, results)

ggplot(df, aes(x=time, y=value, color=k)) + geom_line() + 
  geom_hline(yintercept = colSums(W), linetype='dashed') + 
  facet_grid(.~mstep) + 
  ylab("column norm") + 
  theme_bw() + theme(legend.position = "none")

ggplot(df, aes(x=time, y=value, color=k)) + geom_line() + 
  geom_hline(yintercept = colSums(W), linetype='dashed') + 
  facet_grid(mstep~.) + 
  ylab("column norm") + 
  theme_bw() + theme(legend.position = "none")
