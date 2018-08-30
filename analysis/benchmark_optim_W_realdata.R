# Script used in AISTATS paper / Real dataset
# Demonstrate that W converges faster with updates that only use C (no H)
devtools::load_all()


library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(profvis)

#set.seed(1)
alpha <- 1
beta <- 1
maxiters = 1000
nsamples = 100

# MNIST dataset
# N <- 100
# V <- t(load_mnist())
# F <- dim(V)[1]
# K <- 15
# V <- V[,sample(ncol(V),N)]

# NIPS (full dataset)
# https://archive.ics.uci.edu/ml/datasets/bag+of+words
# https://archive.ics.uci.edu/ml/machine-learning-databases/bag-of-words/
V <- load_nips()
F <- dim(V)[1]
N <- dim(V)[2]
N <- 1000
K <- 10
V <- V[,sample(ncol(V),N)]

#initW <- matrix(runif(F*K), nrow = F, ncol = K)
initW <- replicate(K, rowMeans(V)/K) # Initialize in the center of potential solutions

times <- list()
start <- Sys.time()
res_C <- nmf_mmle_mcem(V, initW, maxiters=maxiters, nsamples, algo='gibbs', mstep="C", alpha=alpha, beta=beta)
w_traces_C <- res_C$W_traces
end <- Sys.time()
times$C <- end-start

res_Csaem <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="C", saem=TRUE, alpha=alpha, beta=beta)
w_traces_Csaem <- res_Csaem$W_traces
times$Csaem <- res_Csaem$times

start <- Sys.time()
res_CH <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="CH", alpha=alpha, beta=beta)
w_traces_CH <- res_CH$W_traces
end <- Sys.time()
times$CH <- end-start

res_CHsaem <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="CH", saem=TRUE, alpha=alpha, beta=beta)
w_traces_CHsaem <- res_CHsaem$W_traces
times$CHsaem <- res_CHsaem$times
save.image("experiment.RData")


# Show evolution of learned dicts
# par(mfrow=c(3,3), pty='s')
# atoms <- 1:9
# for(tr in seq(1,maxiters,by=10)){
#   for (i in atoms){
#     W_sample <- w_traces_C[tr,,atoms]
#     W_sample <- W_sample[,order(-colSums(W_sample))]
#     maxval <- max(W_sample)
#     show_digit(W_sample[,i], zlim=c(0,maxval))
#   }
#   title(tr)
# }

# Show final learned of dict 
# tr <- maxiters
# atoms <- 1:K
# W_sample <- w_traces_C[tr,,atoms]
# W_sample <- W_sample[,order(-colSums(W_sample))]
# maxval <- max(W_sample)
# 
# par(mfrow=c(3,3), pty='s')
# for (i in atoms){
#   show_digit(W_sample[,i], zlim=c(0,maxval))
# }

plot(sort(colSums(res_C$W), decreasing = TRUE)) # plot column norms
plot(sort(colSums(res_CH$W), decreasing = TRUE)) # plot column norms

# Plot traces of column norms (figures used in the AISTATS paper)
################################################################################
margins <- c(1,3)

df_C <- as.data.frame(apply(res_C$W_traces, margins, sum)) %>%
          mutate(mstep = "MCEM-C", iteration=1:maxiters, time=res_C$times) %>%
          gather(k, value, -mstep, -iteration, -time)

df_CH <- as.data.frame(apply(res_CH$W_traces, margins, sum)) %>%
          mutate(mstep = "MCEM-CH", iteration=1:maxiters, time=res_CH$times) %>%
          gather(k, value, -mstep, -iteration, -time)

df_CHsaem <- as.data.frame(apply(res_CHsaem$W_traces, margins, sum)) %>%
          mutate(mstep = "SAEM CH", iteration=1:maxiters, time=res_CHsaem$times) %>%
          gather(k, value, -mstep, -iteration, -time)

df_Csaem <- as.data.frame(apply(res_Csaem$W_traces, margins, sum)) %>%
          mutate(mstep = "SAEM C", iteration=1:maxiters, time=res_Csaem$times) %>%
          gather(k, value, -mstep, -iteration, -time)

df <- rbind(df_C, df_CH)
df <- rbind(df_C, df_CH, df_CHsaem, df_Csaem)


            
ggplot(df, aes(x=iteration, y=value, color=k)) + geom_line() + facet_grid(mstep~.)  + 
  theme_bw() + theme(legend.position = "none", strip.background = element_blank()) + 
  xlab("iteration") + ylab("norm of columns") #+  xlim(0,100)

ggplot(df, aes(x=time, y=value, color=k)) + geom_line() + facet_grid(mstep~.)  + 
  theme_bw() + theme(legend.position = "none", strip.background = element_blank()) + 
  xlab("CPU time (s)") + ylab("norm of columns") #+  xlim(0,100)

# no color
ggplot(df, aes(x=time, y=value, group=k)) + geom_line() + facet_grid(.~mstep)  + 
  theme_bw() + theme(legend.position = "none", strip.background = element_blank()) + 
  xlab("CPU time (s)") + ylab("norm of columns") #+  xlim(0,100)
  ggsave("figreal_nips_2.eps", device = NULL, path = NULL,
         scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
         dpi = 300, limitsize = TRUE)
