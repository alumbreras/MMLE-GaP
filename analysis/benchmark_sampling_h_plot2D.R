# Draws synthetic data, infers the posterior h and plot the traces
# 
library(gtools)
library(gtools)
library(dplyr)
library(tidyr)
library(fields)
library(rstan)
library(coda)

alpha <- 1
beta <- 1
K <- 2
N <- 10
F <- 1000
nsamples <- 500
burnin <- floor(nsamples/2)
# Data -------------------------------------------------------------------------

alphadir <- 100 # high: uniform vectors # low (<1): sparse vectors
#alphadir <- 1 # high: uniform vectors # low (<1): sparse vectors

# inference easier over W with ortogonal rows
W <- rdirichlet(F, alpha = rep(alphadir,K))
W <- W*1

W <- matrix(runif(F*K), ncol=K)
W <- abs(matrix(rnorm(F*K, 1, 10), ncol=K))
k1 <- rep(1:100)
k2 <- rep(100:1)
W <- cbind(k1, k2)
W <- cbind(k1, k1)
W <- W*1

# Generate h from a Gamma and v = Poisson(Wh)
h <- rgamma(K, shape=10, rate=1)
h <- t(t(rep(1,K)))
v <- rpois(F, lambda = W%*%h)



# Plot posterior p(h | ) -------------------------------------------------------
Ngrid <-100
valmax <- 10
h1 <- seq(0.000001,valmax,length.out=Ngrid)
h2 <- seq(0.000001,valmax,length.out=Ngrid)
xygrid <- expand.grid(x=h1, y=h2)
df.xyz <- xygrid %>% 
  rowwise %>% 
  mutate(z = posterior_h(h = c(x,y), v=v, W=W,
                         alpha=alpha, beta=1, log=TRUE)) %>%
  ungroup
z <- df.xyz$z-max(df.xyz$z)
z <- exp(df.xyz$z-max(df.xyz$z))
z <- matrix(z, nrow=sqrt(length(z)))

# Plot on screen
par(pty="s")
plot.new()
#contour(h1, h2, z, add=TRUE, nlevels = 5)
image.plot(z)
contour(h1, h2, z, add = TRUE)

if(FALSE){
  # Plot in file (for the paper)
  par(pty="s")
  plot.new()
  setEPS()
  postscript(paste0("posterior_hard.eps"), width = 9, height = 9)
  image.plot(z)
  contour(h1, h2, z, add = TRUE)
  dev.off()
}
#3D log
z <- df.xyz$z-max(df.xyz$z)
z <- exp(df.xyz$z-max(df.xyz$z))
z <- matrix(z, nrow=sqrt(length(z)))
nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))# Samplers ---------------------------------------------------------------------
zcol  = cut(z, nbcol)
persp3d(h1, h2, z, col=color[zcol])

# initial sample for h
h_init <- rgamma(K, shape=alpha, rate=1) # K-length initialization of h column sampler
all.traces <- list()
times <- list()

# Metropolis -------------------------------------------------------------------
start <- Sys.time()
#traces <- sample_metropolis_gaussian_h(v, W, h_init, alpha=1, beta=1, step=0.01, iter=nsamples)
traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=1, beta=1, step=0.05, iter=nsamples)
cat("\nMH rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$metropolis <- traces
end <- Sys.time()
times$metropolis <- end-start

# Gibbs ----------------------------------------------------------------------
start <- Sys.time()
traces <- sample_gibbs_cpp(v, W, h_init, alpha=1, beta=1, iter=nsamples)$h_n_samples
all.traces$gibbs <- traces
end <- Sys.time()
times$gibbs <- end-start

# HMC --------------------------------------------------------------------------
start <- Sys.time()
traces <- sample_hamiltonian_h(v, W, h_init, epsilon = 0.001, L=10, iter=nsamples)
cat("\nHMC rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$hmc <- traces
end <- Sys.time()
times$hmc <- end-start


# Hamiltonian NUTS STAN---------------------------------------------------------
model <- stan_model("MMLE_H.stan")
start <- Sys.time()
#fit <- rstan::sampling(model, data=list(K=K, F=F, W=W, v=v), iter=nsamples, chains=1)
fit <- rstan::sampling(model, data=list(K=K, F=F, W=W, v=v), iter=nsamples, chains=1, 
                verbose=TRUE, algorithm = "NUTS", init=list(list(h=h_init)))
traces <- as.data.frame(rstan::extract(fit, permuted = FALSE, inc_warmup=TRUE)[,1,1:length(h_init)])
traces <- as.matrix(traces)
cat("\nNUTS rejection (after burnin):", mean(duplicated(traces[-(1:burnin),])))
all.traces$nuts <- traces
end <- Sys.time()
times$nuts <- end-start

# Plot samples over posterior contour ------------------------------------------
xmax <- max(do.call(rbind, all.traces)[,1])
ymax <- max(do.call(rbind, all.traces)[,2])
valmax <- max(do.call(rbind, all.traces))
valmin <- min(do.call(rbind, all.traces))
for(i in 1:length(all.traces)){
  samples <- all.traces[[i]]
  method <- names(all.traces)[i]
  
  xmax <- max(samples[,1])
  ymax <- max(samples[,2])
  xmin <- min(samples[,1])
  ymin <- min(samples[,2])
  valmax <- max(c(xmax, ymax))
  valmin <- min(c(xmin, ymin))
  
  # Plot true likelihood p(V | W)
  Ngrid <-100
  h1 <- seq(0,valmax,length.out=Ngrid)
  h2 <- seq(0,valmax,length.out=Ngrid)
  xygrid <- expand.grid(x=h1, y=h2)
  df.xyz <- xygrid %>% 
    rowwise %>% 
    mutate(z = posterior_h(h = c(x,y), v=v, W=W,
                           alpha=alpha, beta=1, log=TRUE)) %>%
    ungroup
  z <- exp(df.xyz$z-max(df.xyz$z))
  z <- matrix(z, nrow=sqrt(length(z)))
  
  # Plot on screen
  par(pty="s")
  plot.new()
  plot(samples, pch=19, cex=0.2, col='red', xlab="", ylab="")
  points(t(h), col="blue", pch=15, cex=1)
  contour(h1, h2, z, add=TRUE, nlevels = 20)
  title(method)

  if(FALSE){
    # Plot in file (for the paper)
    par(pty="s")
    plot.new()
    setEPS()
    postscript(paste0("samples_2d_", names(all.traces)[i],".eps"), width = 9, height = 9)
    plot(samples, pch=19, cex=0.2, col='red', xlab="", ylab="")
    contour(h1, h2, z, add=TRUE, nlevels = 20)
    dev.off()
  }
}


# Plots traces and autocorrelations --------------------------------------------
mytheme <-  theme_classic() + 
            theme(legend.position="none", axis.text = element_text(size=10),
                  panel.border = element_rect(colour = "black", fill=NA, size=0.5))
            

mywidth <- 10
myheight <- 10
all.ess <- list()
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
    geom_hline(data=df.ground_truth, aes(yintercept=h)) + 
    ylab("") +
    mytheme + ggtitle(method)
  print(p)
  ggsave(paste0("traces_2d_", names(all.traces)[i],".eps"), device="eps", 
         width=mywidth, height=myheight, units="cm")
  
  # Plot autocorrelations
  traces <- all.traces[[i]]
  mcmc.traces <- mcmc(traces)
  autocorr.plot(mcmc.traces, lag.max = 100, auto.layout = TRUE)
  ess <- effectiveSize(mcmc.traces)
  title(method)
  print(ess)

  acf_h1 <- acf(mcmc.traces[,1])$acf[,,1]
  acf_h2 <- acf(mcmc.traces[,2])$acf[,,1]
  df <- data.frame(lag = 1:length(acf_h1), h1=acf_h1, h2=acf_h2)
  df <- gather(df, h, ACF, h1, h2)
  p <- ggplot(df, aes(x=lag, y=ACF, group=h, color=h))  + 
    geom_point() + 
    geom_segment(mapping = aes(xend = lag, yend = 0))+
    geom_hline(aes(yintercept = 0))+
    ylab("") +
    mytheme + ggtitle(method)
  print(p)
  ggsave(paste0("ACF_2d_", names(all.traces)[i],".eps"), device="eps", 
         width=mywidth, height=myheight, units="cm")  
  all.ess[[method]] <- mean(ess)
  
  # Plot trace likelihood
  likelihoods <- apply(all.traces[[method]], 1, function(x) posterior_h(x, v, W, alpha, beta, log=TRUE))
  df.likelihoods <- data.frame(sample=1:length(likelihoods), likelihood = likelihoods )
  p <- ggplot(df.likelihoods, aes(x=sample, y=likelihood)) + 
    geom_line() + mytheme  + ggtitle(method)
  print(p)
  ggsave(paste0("likelihood_2d_", names(all.traces)[i],".eps"), device="eps", 
         width=mywidth, height=myheight, units="cm") 
}

res <- list()
for(method in names(all.traces)){
  likelihoods <- apply(all.traces[[method]], 1, function(x) posterior_h(x, v, W, alpha, beta, log=TRUE))
  res[[method]] <- data.frame(method = method, sample=1:length(likelihoods), likelihood = likelihoods )
}
df.likelihoods <- do.call("rbind", res)
p <- ggplot(df.likelihoods, aes(x=sample, y=likelihood, group=method, color=method)) + 
  geom_line() + theme_bw()
print(p)