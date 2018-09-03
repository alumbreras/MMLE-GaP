# Demonstrates that W converges faster with updates that only use C (no H)
devtools::load_all()
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(latex2exp)

# generate synthetic data from a GaP model ####
maxiters = 100
nsamples = 250

# ICML 
F <- 4
N <- 100
K_star <- 2
K <- 3
alpha_star <- 1
beta_star <- 1
alpha <- 1
beta <- 1

W <- t(rdirichlet(K_star, alpha = rep(1,F))) # ICML
H <- matrix(rgamma(K_star*N, shape=alpha_star, rate=beta_star), ncol=N)
V <- matrix(NA, nrow=F, ncol=N)
for(f in 1:F){
  for(n in 1:N){
    V[f,n] <- rpois(1, lambda = W[f,] %*% H[,n])
  }
}

cat("Mean/Var =", mean(c(V))/var(c(V)))

# Initialize in the center of potential solutions
initW <- replicate(K, rowMeans(V)/K) 

# Algorithms -------------------------------------------------------------------
#V <- V[,colSums(V)>0] # remove empty columns
times <- list()
res_C <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="C", alpha=alpha, beta=beta)
w_traces_C <- res_C$W_traces
times$C <- res_C$times

res_CH <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="CH", alpha=alpha, beta=beta)
w_traces_CH <- res_CH$W_traces
times$CH <- res_CH$times

nsamples <- 20
res_Csaem <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="C", saem=TRUE, alpha=alpha, beta=beta)
w_traces_Csaem <- res_Csaem$W_traces
times$Csaem <- res_Csaem$times

res_CHsaem <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="CH", saem=TRUE, alpha=alpha, beta=beta)
w_traces_CHsaem <- res_CHsaem$W_traces
times$CHsaem <- res_CHsaem$times

res_H <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="H", alpha=alpha, beta=beta)
w_traces_H <- res_H$W_traces
times$H <- res_H$times

res_Hu <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="H", mult_updates=100, alpha=alpha, beta=beta)
w_traces_Hu <- res_Hu$W_traces
times$Hu <- res_Hu$times

res_Z <- nmf_mmle_mcem(V[,colSums(V)>0], initW, maxiters, nsamples, algo='gibbs_z', alpha=alpha, beta=beta)
w_traces_Z <- res_Z$W_traces
times$Z <- res_Z$times

plot(sort(colSums(res_C$W), decreasing = TRUE)) # plot column norms
plot(sort(colSums(res_CH$W), decreasing = TRUE)) # plot column norms
plot(sort(colSums(res_Z$W), decreasing = TRUE)) # plot column norms

# Plot evolution of column norms (paper)
################################################################################
margins <- c(1,3) # for column norms, sum over 2nd dimension (K)
# margins <- c(1,2) # for row norms, sum over 3nd dimension (F)

df.colsums_CH <- as.data.frame(apply(res_CH$W_traces, margins, sum)) %>%
                  mutate(mstep = "MCEM CH", iteration=1:maxiters, time=res_CH$times) %>%
                  gather(k, value, -mstep, -iteration, -time)

df.colsums_C <- as.data.frame(apply(res_C$W_traces, margins, sum)) %>%
                  mutate(mstep = "MCEM C", iteration=1:maxiters, time=res_C$times) %>%
                  gather(k, value, -mstep, -iteration, -time)

df.colsums_CHsaem <- as.data.frame(apply(res_CHsaem$W_traces, margins, sum)) %>%
                  mutate(mstep = "SAEM CH", iteration=1:maxiters, time=res_CHsaem$times) %>%
                  gather(k, value, -mstep, -iteration, -time)

df.colsums_Csaem <- as.data.frame(apply(res_Csaem$W_traces, margins, sum)) %>%
                  mutate(mstep = "SAEM C", iteration=1:maxiters, time=res_Csaem$times) %>%
                  gather(k, value, -mstep, -iteration, -time)

df.colsums_H <- as.data.frame(apply(res_H$W_traces, margins, sum)) %>%
                  mutate(mstep = "MCEM H", iteration=1:maxiters, time=res_H$times) %>%
                  gather(k, value, -mstep, -iteration, -time)

df.colsums_Hu <- as.data.frame(apply(res_Hu$W_traces, margins, sum)) %>%
                  mutate(mstep = "MCEM Hu", iteration=1:maxiters, time=res_Hu$times) %>%
                  gather(k, value, -mstep, -iteration, -time)

df.colsums_Z <- as.data.frame(apply(res_Z$W_traces, margins, sum)) %>%
                mutate(mstep = "MCEM Z", iteration=1:maxiters, time=res_Z$times) %>%
                gather(k, value, -mstep, -iteration, -time)

df <- rbind(df.colsums_CH, df.colsums_C)
df <- rbind(df.colsums_CH, df.colsums_C, df.colsums_H)
df <- rbind(df.colsums_CH, df.colsums_C, df.colsums_H, df.colsums_Hu)
df <- rbind(df.colsums_CH, df.colsums_C, df.colsums_H, df.colsums_Z)
df <- rbind(df.colsums_CH, df.colsums_C, df.colsums_CHsaem, df.colsums_Csaem)
  
text <- TeX(sprintf('$F=%d, K=%d, K*=%d, N=%d, \\alpha = %.1f, \\beta = %.1f, samples/iter=%d, iter=%d$', 
                    F, K, K_star, N, alpha, beta, nsamples, maxiters))

ggplot(df, aes(x=iteration, y=value, group=k)) + 
  geom_line() + 
  facet_grid(mstep~.) + 
  geom_hline(yintercept = as.vector(colSums(W)), linetype='dashed') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) + 
  xlab("iteration") + ylab("norm of columns") +
  ggtitle(text)

ggplot(df, aes(x=time, y=value, group=k)) + 
  geom_line() + 
  facet_grid(mstep~.) + 
  geom_hline(yintercept = as.vector(colSums(W)), linetype='dashed') +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white")) + 
  xlab("CPU time (s)") + ylab("norm of columns") +
  ggtitle(text)

ggsave("figtoy.png", device = NULL, path = NULL,
              scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
              dpi = 300, limitsize = TRUE)
       
#ggsave("figtoy_realistic.eps", device = NULL, path = NULL,
#       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
#       dpi = 300, limitsize = TRUE)

###########################################################
# Plot every trace W_ij trace (for small F, K)
###########################################################
res_CH$W_traces2 <- res_CH$W_traces
res_CH$W_traces2[,,1] <- res_CH$W_traces[,,2]
res_CH$W_traces2[,,2] <- res_CH$W_traces[,,1]


df.traces_C <- as.data.frame(res_C$W_traces) %>%
  mutate(mstep = "MCEM C", iteration=1:maxiters, time=res_C$times) %>%
  gather(variable, w, -mstep, -iteration, -time)

df.traces_CH <- as.data.frame(res_CH$W_traces2) %>%
  mutate(mstep = "MCEM CH", iteration=1:maxiters, time=res_CH$times) %>%
  gather(variable, w, -mstep, -iteration, -time)

df.traces_H <- as.data.frame(res_H$W_traces) %>%
  mutate(mstep = "MCEM H", iteration=1:maxiters, time=res_H$times) %>%
  gather(variable, w, -mstep, -iteration, -time)

df.traces_Z <- as.data.frame(res_Z$W_traces) %>%
  mutate(mstep = "MCEM Z", iteration=1:maxiters, time=res_Z$times) %>%
  gather(variable, w, -mstep, -iteration, -time)

df <- rbind(df.traces_CH, df.traces_C)
df <- rbind(df.traces_CH, df.traces_C, df.traces_H)
df <- rbind(df.traces_CH, df.traces_C, df.traces_H, df.traces_Z)

# iterations
ggplot(df, aes(x=iteration, y=w, color=variable)) + geom_line() + 
  facet_grid(mstep~.) + geom_hline(yintercept = as.vector(W), linetype='dashed'  ) + 
  theme_bw() + theme(legend.position = "none")
# time
ggplot(df, aes(x=time, y=w, color=variable)) + geom_line() + 
  facet_grid(mstep~.) + geom_hline(yintercept = as.vector(W), linetype='dashed') + 
  theme_bw() + theme(legend.position = "none") + xlab("CPU time (s)")

# facet_grid side by side
ggplot(df, aes(x=time, y=w, color=variable)) + geom_line() + 
  facet_grid(.~mstep) + geom_hline(yintercept = as.vector(W), linetype='dashed') + 
  theme_bw() + theme(legend.position = "none") + xlab("CPU time (s)")

# Color by columns (F=4, K=3)
df$column <- NA
df$column[df$variable == "V1" | df$variable == "V2" | df$variable == "V3" | df$variable == "V4"] <- "k1"
df$column[df$variable == "V5" | df$variable == "V6" | df$variable == "V7" | df$variable == "V8"] <- "k2"
df$column[df$variable == "V9" | df$variable == "V10" | df$variable == "V11" | df$variable == "V12"] <- "k3"
ggplot(df, aes(x=time, y=w, color=column, group=variable)) + geom_line() + 
  facet_grid(mstep~.) + geom_hline(yintercept = as.vector(W), linetype='dashed') + 
  theme_bw() + theme(legend.position = "none") + xlab("CPU time (s)")
ggsave("figtoy6_traces.png", device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)


# Likelihoods
#############################################################
#likelihood_V_W_cpp(V,  W, alpha=alpha, beta=1)
#likelihood_V_W(V,  W, alpha=alpha, beta=1, log=TRUE, vectorizable=TRUE) # test
df.likelihood_CH <- data.frame(time = times$CH, iteration = 1:maxiters, mstep = "CH")
df.likelihood_CH$z <- apply(w_traces_CH, 1, function(vars) {
                        likelihood_V_W_cpp(V,  W = matrix(vars, nrow=F, ncol=K), 
                                       alpha=alpha, beta=alpha)})

df.likelihood_C <- data.frame(time = times$C, iteration = 1:maxiters, mstep = "C")
df.likelihood_C$z <- apply(w_traces_C, 1, function(vars) {
                        likelihood_V_W_cpp(V,  W = matrix(vars, nrow=F, ncol=K), 
                                             alpha=alpha, beta=alpha)})

df.likelihood_H <- data.frame(time = times$H, iteration = 1:maxiters, mstep = "H")
df.likelihood_H$z <- apply(w_traces_H, 1, function(vars) {
  likelihood_V_W_cpp(V,  W = matrix(vars, nrow=F, ncol=K), 
                     alpha=alpha, beta=alpha)})

df.likelihood_Z <- data.frame(time = times$H, iteration = 1:maxiters, mstep = "Z")
df.likelihood_Z$z <- apply(w_traces_Z, 1, function(vars) {
  likelihood_V_W_cpp(V,  W = matrix(vars, nrow=F, ncol=K), 
                     alpha=alpha, beta=alpha)})

df <- rbind(df.likelihood_CH, df.likelihood_C)
df <- rbind(df.likelihood_CH, df.likelihood_C, df.likelihood_H)
df <- rbind(df.likelihood_CH, df.likelihood_C, df.likelihood_H, df.likelihood_Z)

df <- select(df, iteration, time, mstep, z)
df$mstep <- as.factor(df$mstep)
ggplot(df, aes(x=iteration, y=z, group=mstep, color=mstep)) + geom_line() + ylab("loglikelihood")
ggplot(df, aes(x=time, y=z, group=mstep, color=mstep)) + geom_line() + ylab("loglikelihood")
ggsave("figtoy6_likelihoods.png", device = NULL, path = NULL,
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)