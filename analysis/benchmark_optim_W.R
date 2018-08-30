# Demonstrate theat W converges faster with updates that only use C (no H)


library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(cowplot)

set.seed(1)

maxiters = 400
nsamples = 250

# generate synthetic data from a GaP model ####
F <- 1
K <- 2
#F <- 100
#K <- 10
N <- 100
W <- matrix(c(1,0), nrow=F,  ncol=K) # F x K matrix
F <- nrow(W)
K <- ncol(W)
H <- matrix(rgamma(K*N, shape=alpha, rate=1), ncol=N) # K x N matrix
V <- matrix(NA, nrow=F, ncol=N)
for(f in 1:F){
  for(n in 1:N){
    V[f,n] <- rpois(1, lambda = W[f,] %*% H[,n])
  }
}

initW <- matrix(1, nrow = F, ncol = K) # K initially overestimated

res <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="C")
w_traces_C <- res$W_traces[,1,]

res <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="H")
w_traces_H <- res$W_traces[,1,]

res <- nmf_mmle_mcem(V, initW, maxiters, nsamples, algo='gibbs', mstep="CH")
w_traces_CH <- res$W_traces[,1,]

plot(sort(colSums(res$W), decreasing = TRUE)) # plot column norms

df <- as.data.frame(w_traces_H)
df$iteration <- 1:nrow(df)
df$sum <- df$V1 + df$V2
df$mstep <- "H"
dfH <- tidyr::gather(df, variable, value, -iteration, -mstep)
#dfH <- df[complete.cases(dfH),]

df <- as.data.frame(w_traces_C)
df$iteration <- 1:nrow(df)
df$sum <- df$V1 + df$V2
df$mstep <- "C"
dfC <- tidyr::gather(df, variable, value, -iteration, -mstep)

df <- as.data.frame(w_traces_CH)
df$iteration <- 1:nrow(df)
df$sum <- df$V1 + df$V2
df$mstep <- "CH"
dfCH <- tidyr::gather(df, variable, value, -iteration, -mstep)

df <- rbind(dfH, dfC, dfCH)
df$mstep <- as.factor(df$mstep)
ggplot(df, aes(x=iteration, y=value, color=variable)) + geom_line() + facet_wrap(~mstep)
ggplot(df, aes(x=iteration, y=value, color=variable)) + geom_line() + facet_grid(mstep~.)

# Louis' plot
if(K==2){
  N <- 50
  w11 <- seq(0,10,length.out=N)
  w12 <- seq(0,10,length.out=N)
  xygrid <- expand.grid(x=w11, y=w12)
  
  df.xyz <- xygrid %>% 
    rowwise %>% 
    mutate(z = likelihood_V_W(V,  W = matrix(c(x,y), nrow=1, ncol=2), 
                              alpha=alpha, beta=1, log= FALSE, vectorizable=TRUE)) %>%
    ungroup
  
  df.xyz %>% arrange(-z)
  MLE_probs <- df.xyz[which.max(df.xyz$z),]
  
  z <- df.xyz$z
  z <- matrix(z, nrow=sqrt(length(z)))
  image.plot(z)
  contour(z, add = TRUE)
  
  par(pty="s", mfrow=c(1,2))
  for(i in 1:10){
    plot(w_traces_C[i,1], w_traces_C[i,2], pch=19, cex=0.2, col='red', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
    points(w_traces_H[i,1], w_traces_H[i,2], pch=19, cex=0.2, col='blue', xlab="", ylab="")
    title(i)
  }
    
  #################""
  par(pty="s", mfrow=c(1,3))
  plot(w_traces_H, pch=19, cex=0.2, col='red', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  lines(w_traces_C, pch=19, cex=0.2, col='blue', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  lines(w_traces_CH, pch=19, cex=0.2, col='green', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  
  par(pty="s", mfrow=c(1,3))
  plot(w_traces_C, pch=19, cex=0.2, col='red', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  plot(w_traces_H, pch=19, cex=0.2, col='red', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  plot(w_traces_CH, pch=19, cex=0.2, col='red', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  
  for(from in seq(1, maxiters-10, by=10)){
  to <- from+10
  plot(w_traces_H[from:to,], pch=19, cex=0.2, col='red', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  lines(w_traces_C[from:to,], pch=19, cex=0.2, col='blue', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  lines(w_traces_CH[from:to,], pch=19, cex=0.2, col='green', xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), type='o')
  readline()
  }

}


par(mfrow=c(1,1))
error_H <- apply(w_traces_H, 1, function(x) sum((x-W)^2))
error_C <- apply(w_traces_C, 1, function(x) sum((x-W)^2))
error_CH <- apply(w_traces_CH, 1, function(x) sum((x-W)^2))
maxerror <- max(error_H, error_C, error_CH)
plot(error_H, pch=19, cex=0.2, col='red', xlab="", ylab="", type='o', ylim=c(0,maxerror))
lines(error_C, pch=19, cex=0.2, col='blue', xlab="", ylab="", type='o', ylim=c(0,maxerror))
lines(error_CH, pch=19, cex=0.2, col='green', xlab="", ylab="", type='o', ylim=c(0,maxerror))


df.likelihood_W_H <- as.data.frame(w_traces_H) %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(V1,V2), nrow=1, ncol=2), 
                            alpha=alpha, beta=1, log= TRUE, vectorizable=TRUE)) %>%
  ungroup

df.likelihood_W_C <- as.data.frame(w_traces_C) %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(V1,V2), nrow=1, ncol=2), 
                            alpha=alpha, beta=1, log= TRUE, vectorizable=TRUE)) %>%
  ungroup

df.likelihood_W_CH <- as.data.frame(w_traces_CH) %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(V1,V2), nrow=1, ncol=2), 
                            alpha=alpha, beta=1, log= TRUE, vectorizable=TRUE)) %>%
  ungroup
plot(df.likelihood_W_H$z, type="l", col="red")
lines(df.likelihood_W_C$z, type="l", col="blue")
lines(df.likelihood_W_CH$z, type="l", col="green")