# Benchmark multiple Monte Carlo samplers for a column pf H 
# using all the dictionary and the corresponding column of V.
# Do it using a grid of parameters.

# Plot, RMSE (after burnin)
# Plot effective sample size

library(dplyr)
library(tidyr)
library(fields)
library(rgl)
library(gtools)
library(rstan)
library(coda)
library(profvis)
library(mcmcse)
library("pixmap")

plot_vector_image <- function(vimage){
  side <- sqrt(length(vimage))
  mimage <- matrix(vimage, ncol=side)
  image(t(mimage)[,side:1], col = gray.colors(256))
}

nsamples <- 500
burnin <- floor(nsamples/2)

df.results <-  rep(NA, 8)
df.results <-  data.frame(alphadir=NA, alpha=NA, F=NA, K=NA, metropolis=NA, gibbs=NA, hamiltonian=NA, nuts=NA)
df.results.ess <- df.results
df.results.time <- df.results

alphadir <- 1
df.summary.last <- data.frame()


mhstep1 <- 0.01
mhstep2 <- 0.05
mhstep3 <- 0.1
hmc1 <- 15
hmc2 <- 25

# Original image
h <- read.pnm(file = "../data/cameraman/cameraman64.pgm")
h <- as.vector(h@grey*255)


# Load blurring matrix
W <- read.csv("../data/cameraman/W_blurring.csv", header = FALSE)
W <- t(as.matrix(W))
F <- dim(W)[1]
K <- dim(W)[2]

Wh <- W%*%h
v <- rpois(F, lambda = W%*%h)
par(mfrow=c(1,3))
plot_vector_image(h)
plot_vector_image(Wh)
plot_vector_image(v)

# initial sample for h
h_init <- rgamma(K, shape=10, rate=beta) # K-length initialization of h column sampler
# ******************************************************************************
all.traces <- list()
all.likes <- list()
times <- list()

# Variational
# use this values to create an h_init (the MAP of this variational distribution)
alpha_var <- runif(K)
beta_var <- runif(K)
res_vb <- optim_variational(v, W, maxiters = 100, alpha=alpha, beta=beta)

res_vb <- variational_ch_cpp(v, W, alpha_var=alpha_var, beta_var=beta_var, 
                            alpha=alpha, beta=beta, maxiters = 1000)
plot(res_vb$lp, main="lower bound")
alpha_var <- res_vb$q_h_alpha
beta_var <- res_vb$q_h_beta
h_hat <- c(alpha_var/beta_var) # Expectation
plot_vector_image(h_hat)
alpha_prior <- 0.001
beta_prior <- 0.001
plot(res_vb$lower_bound)
plot(exp(res_vb$lower_bound-max(res_vb$lower_bound)))

# Metropolis -------------------------------------------------------------------
start <- Sys.time()
traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=alpha_prior, beta=beta_prior, step=mhstep1, iter=nsamples)
end <- Sys.time()
cat("\nMetropolis rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$metropolis1 <- traces
all.likes$metropolis1 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
times$metropolis1 <- end-start
cat("\nMetropolis took", end-start, "secs")

start <- Sys.time()
traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=alpha_prior, beta=beta_prior, step=mhstep2, iter=nsamples)
end <- Sys.time()
cat("\nMetropolis rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$metropolis2 <- traces
all.likes$metropolis2 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
times$metropolis2 <- end-start
cat("\nMetropolis took", end-start, "secs")

start <- Sys.time()
traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=alpha_prior, beta=beta_prior, step=mhstep3, iter=nsamples)
end <- Sys.time()
cat("\nMetropolis rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$metropolis3 <- traces
all.likes$metropolis3 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
times$metropolis3 <- end-start
cat("\nMetropolis took", end-start, "secs")

burnin <- floor(dim(traces)[1]/2)
h_hat <- colMeans(traces[-(1:burnin),])
h_hat <- traces[nsamples,]
plot_vector_image(h_hat)


h_hat <- traces[nsamples,]
plot_vector_image(h_hat)

# Gibbs ------------------------------------------------------------------------
start <- Sys.time()
traces <- sample_gibbs_cpp_scalable(v, W, h_init, alpha=alpha_prior, beta=beta_prior, iter=nsamples)$h_n_samples
end <- Sys.time()
all.traces$gibbs <- traces
all.likes$gibbs <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
times$gibbs <- end-start
cat("\nGibbs took", end-start, "secs")

burnin <- floor(dim(traces)[1]/2)
h_hat <- colMeans(traces[-(1:burnin),])
h_hat <- traces[nsamples,]
plot_vector_image(h_hat)

# Hamiltonian -----------------------------------------------------------------
start <- Sys.time()
traces <- sample_hmc_cpp(v, W, h_init, alpha = alpha_prior, beta = beta_prior, epsilon = 0.001, L=hmc1, iter=nsamples)
end <- Sys.time()
cat("\nHamiltonian rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$hmc1 <- traces
all.likes$hmc1 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
times$hmc1 <- end-start
cat("\nHMC took", end-start, "secs")

# Hamiltonian -----------------------------------------------------------------
start <- Sys.time()
traces <- sample_hmc_cpp(v, W, h_init, alpha = alpha_prior, beta = beta_prior, epsilon = 0.001, L=hmc2, iter=nsamples)
end <- Sys.time()
cat("\nHamiltonian rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$hmc2 <- traces
all.likes$hmc2 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
times$hmc2 <- end-start
cat("\nHMC took", end-start, "secs")

burnin <- floor(dim(traces)[1]/2)
h_hat <- colMeans(traces[-(1:burnin),])
h_hat <- traces[nsamples,]
plot_vector_image(h_hat)

# NUTS (Stan)---------------------------------------------------------
model <- rstan::stan_model("MMLE_H2.stan")
#stanc("MMLE_H.stan")
#cat(stanc("MMLE_H.stan")$cppcode)
start <- Sys.time()
fit <- rstan::sampling(model, data=list(K=K, F=F, W=W, v=v, alpha=alpha_prior, beta=beta_prior), 
                       iter=nsamples, chains=1, 
                       verbose=TRUE, algorithm = "NUTS", init=list(list(h=h_init)))
end <- Sys.time()
traces <- as.data.frame(rstan::extract(fit, permuted = FALSE, inc_warmup=TRUE)[,1,1:length(h_init)])
traces <- as.matrix(traces)
cat("\nNUTS rejection:",  mean(duplicated(traces[-(1:burnin),])))
all.traces$nuts <- traces
all.likes$nuts <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
times$nuts <- end-start
cat("\nNUTS (Stan) took", end-start, "secs")

burnin <- floor(dim(traces)[1]/2)
h_hat <- colMeans(traces[-(1:burnin),])
h_hat <- traces[nsamples,]
plot_vector_image(h_hat)

# Plot images of posterior mean -----------------
par(pty="s", mfrow=c(3,3))
for(i in 1:length(all.traces)){
  method <- names(all.traces)[i]
  traces <- all.traces[[i]][-(1:burnin),]
  plot_vector_image(colMeans(traces))
  title(method)
}


# Store obtained metrics -------------------
for(i in 1:length(all.traces)){
  method <- names(all.traces)[i]
  traces <- all.traces[[i]][-(1:burnin),]
  df <- data.frame(alphadir=alphadir, alpha=alpha, F=F, K=K,
                   method  = method,
                   worst_z = max(abs(geweke.diag(mcmc(traces))$z)),
                   time    = as.numeric(times[[method]], units="secs"),
                   ESS     = mean(effectiveSize(mcmc(traces))),
                   minESS  = min(effectiveSize(mcmc(traces))),
                   RMSE    = sqrt(sum(apply(traces, 1, function(x) t(x-h)%*%(x-h)))/nsamples),
                   maxlike = max(apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha_prior, beta_prior))),
                   accepted= 1 - mean(duplicated(traces))) 
  
  df.summary.last <- rbind(df.summary.last, df)
}


# Iters vs SSE and likelihood
df.traces_quality <- data.frame()
for(i in 1:length(all.traces)){
  method <- names(all.traces)[i]
  traces <- all.traces[[i]]
  nsamples_ <- nrow(traces)
  milsecs <- seq(0,as.numeric(times[[method]], units="secs"), length.out=nsamples_)*1000
  df.traces_quality <- rbind(df.traces_quality, 
                     data.frame(
                        SE = apply(traces, 1, function(x) t(x-h)%*%(x-h)),
                        #Frob = apply(traces, 1, function(x) norm(x-h, type="F")),
                        likelihoodCpp_ = apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha_prior, beta_prior)),
                        likelihoodCpp = apply(traces, 1, function(x) posterior_h_cpp_3(x, v, W, alpha_prior+0.5, beta_prior+0.5)),
                        likelihood = apply(traces, 1, function(x) posterior_h(x, v, W, alpha_prior, beta_prior,log=TRUE)),
                        method     = method, 
                        iteration  = 1:nsamples_, 
                        time=milsecs))
}

df.traces_quality$method <- as.character(df.traces_quality$method)
df.traces_quality[df.traces_quality$method == "metropolis1",]$method <- "MH-1"
df.traces_quality[df.traces_quality$method == "metropolis2",]$method <- "MH-2"
df.traces_quality[df.traces_quality$method == "metropolis3",]$method <- "MH-3"
df.traces_quality[df.traces_quality$method == "gibbs",]$method <- "Gibbs"
df.traces_quality[df.traces_quality$method == "hmc1",]$method <- "HMC-1"
df.traces_quality[df.traces_quality$method == "hmc2",]$method <- "HMC-2"
df.traces_quality[df.traces_quality$method == "nuts",]$method <- "NUTS"

df.traces_quality$method <- as.factor(df.traces_quality$method)

levels <- c("MH-1","MH-2","MH-3","HMC-1","HMC-2","Gibbs","NUTS")
df.traces_quality$method <- factor(df.traces_quality$method, levels = levels)

mytheme <- theme_bw() + theme(aspect.ratio=1/2.5, legend.title = element_blank(), legend.position = "none",
                              strip.background = element_blank(),
                              text = element_text(size=20))

ggplot(df.traces_quality,
       aes(x=iteration, y=likelihoodCpp_, color=method, group=method)) + 
      geom_line() +
      facet_grid(method~.)

  ggplot(df.traces_quality,
         aes(x=iteration, y=SE, color=method, group=method)) + 
  geom_line() + 
  scale_y_log10() +
  facet_grid(method~.)


#####################################################"

# Store all traces ready to be plotted
df.traces <- data.frame()
for(i in 1:length(all.traces)){
  method <- names(all.traces)[i]
  traces <- all.traces[[i]][1:10,]
  nsamples_ <- nrow(traces)
  milsecs <- seq(0,as.numeric(times[[method]], units="secs"), length.out=nsamples_)*1000
  df.traces <- rbind(df.traces, 
                     as.data.frame(traces) %>%
                     mutate(method     = method, 
                            iteration  = 1:nsamples_, 
                            time=milsecs) %>%
                     gather(k, value, -method, -iteration, -time))
}
p <- ggplot(df.traces, aes(x=iteration, y=value, color=k)) + geom_line() + 
  geom_hline(yintercept = h, linetype='dashed') +
  facet_grid(method~.) +
  scale_y_log10()+
  #ylim(0,100) +
  theme_bw() + theme(legend.position = "none")
print(p)


# Store Likelihoods (not used) -------------------------------------------------
res <- list()
for(method in names(all.traces)){
  likelihoods <- apply(all.traces[[method]], 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
  res[[method]] <- data.frame(method = method, sample=1:length(likelihoods), likelihood = likelihoods )
}
###############################################################################
cat("\n\nEnd of experiments. Plotting the results")

df.summary <- df.summary.last
df.summary.backup <- df.summary #just in case
df.summary$time <- as.numeric(df.summary$time)
df.summary$method <- as.character(df.summary$method)

df.summary[df.summary$method == "metropolis1",]$method <- "MH-1"
df.summary[df.summary$method == "metropolis2",]$method <- "MH-2"
df.summary[df.summary$method == "metropolis3",]$method <- "MH-3"
df.summary[df.summary$method == "gibbs",]$method <- "Gibbs"
df.summary[df.summary$method == "hmc1",]$method <- "HMC-1"
df.summary[df.summary$method == "hmc2",]$method <- "HMC-2"
df.summary[df.summary$method == "nuts",]$method <- "NUTS"

df.summary$method <- as.factor(df.summary$method)

levels <- c("MH-1","MH-2","MH-3","HMC-1","HMC-2","Gibbs","NUTS")
df.summary$method <- factor(df.summary$method, levels = levels)

mytheme <- theme_bw() + theme(aspect.ratio=1/2.5, legend.title = element_blank(), legend.position = "none",
                              strip.background = element_blank(),
                              text = element_text(size=20))

mytheme_legend <- theme_bw() + theme(aspect.ratio=1/2.5, legend.title = element_blank(),
                              strip.background = element_blank(),
                              text = element_text(size=20))
# plot.margin=grid::unit(c(0,0,0,0), "mm"))

# RMSE #######
# w.r.t F
ggplot(df.summary,
       aes(x=1, y=ESS, shape=method)) +  
  scale_y_log10()+
  scale_shape_manual(values=1:nlevels(df.summary$method)) +
  geom_point(aes(group=method, shape=method)) + 
  mytheme_legend

ggsave("cameraman_RMSE.eps",
       scale = 1, units = "cm", dpi = 300)

# worst_z #######
ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=worst_z, fill=method, color=F)) +  
  facet_grid(K~.) +
  scale_shape_manual(values=1:nlevels(df.summary$method)) +
  scale_color_manual(values=1:nlevels(df.summary$method), guide=FALSE) +
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) + 
  ylab("maximum Geweke z-score ") + ylim(0,50)+
  mytheme

ggsave("cameraman_worst_Z.eps",
       scale = 1, units = "cm", dpi = 300)


# ESS/time ########
ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=minESS/time, fill=method, color=F)) +  
  facet_grid(K~.) +
  scale_shape_manual(values=1:nlevels(df.summary$method)) +
  scale_color_manual(values=1:nlevels(df.summary$method), guide=FALSE) +
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) +
  scale_y_log10() +
  ylab("minimum ESS/time") +
  mytheme

ggsave("cameraman_esstime.eps",
       scale = 1, units = "cm", dpi = 300)

#*******************************************************************************
# NOT USED IN THE PAPER
#*******************************************************************************

# ESS
############
# w.r.t F
ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=minESS, fill=method)) +  
  facet_grid(K~.) +
  scale_y_log10()+
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) + 
  theme_bw()