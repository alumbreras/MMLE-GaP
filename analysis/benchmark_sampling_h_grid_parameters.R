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

nsamples <- 1000
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

F <- 1000 # for quick tests
K <- 5

# repeat 10 times 
for(xp in 1){
  for(alphadir in c(1)){
    #for(F in c(100, 1000, 10000, 50000)){
    #for(F in c(100, 1000, 10000)){
    #  for(K in c(10,100)){
    for(F in c(1000)){
      for(K in c(10)){
        K = 5
        F = 100
        cat("\n ************************************************************************\n")
        cat(paste("xp:", xp, "alphadir:", alphadir, "F:", F, "K:", K))
        # ******************************************************************************
        # alphadir high: uniform vectors # low (<1): sparse vectors
        # F high: easier inference
        # K high: no effect? (just slower)
        # Data -------------------------------------------------------------------------
        # ******************************************************************************
        case <- "hard"
        if(case == "hard"){
          # Generate H from a Gamma and V = Poisson(WH)
          alpha <- 1
          beta <- 1
          alphadir <- 100
          h <- rgamma(K, shape=1, rate=0.01)
          h <- rgamma(K, shape=10, rate=1)
          h <- seq(1,K*100, length.out = K)
          #h <- 10^(1:K)   #seq(1,K*100, length.out = K)
          W <- rdirichlet(F, alpha = rep(alphadir,K))
          v <- rpois(F, lambda = W%*%h)
        }
        
        # initial sample for h
        h_init <- rgamma(K, shape=10, rate=beta) # K-length initialization of h column sampler
        h_init <- h
        # ******************************************************************************
        all.traces <- list()
        all.likes <- list()
        times <- list()
        
        # VB -------------------------------------------------------------------
        alpha_var <- runif(K)
        beta_var <- runif(K)
        start <- Sys.time()
        res_vb <- variational_ch_cpp(v, W, alpha_var=alpha_var, beta_var=beta_var, 
                                     alpha=alpha, beta=beta, maxiters = 100)
        end <- Sys.time()
        
        var_alpha <- res_vb$q_h_alpha
        var_beta <- res_vb$q_h_beta
        traces <- matrix(NA, nrow = nsamples, ncol = K)
        for(kk in 1:K){
            traces[,kk] <- rgamma(nsamples, var_alpha[kk], rate = var_beta[kk])
        }
        all.traces$VB <- traces
        all.likes$VB <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        times$VB <- end-start
        ellapsed <- as.numeric(end-start, units="secs")
        cat("\nVB took", ellapsed, "secs")
        plot(res_vb$lp)
        plot_traces(traces, h, title="MH")
        
        h_init <- c(var_alpha / var_beta) # Expectation
        
        # Metropolis -------------------------------------------------------------------
        start <- Sys.time()
        traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=1, beta=1, step=0.01, iter=nsamples)
        end <- Sys.time()
        cat("\nMetropolis rejection:",  mean(duplicated(traces[-(1:burnin),])))
        all.traces$MH <- traces
        all.likes$MH <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        times$MH <- end-start
        ellapsed <- as.numeric(end-start, units="secs")
        cat("\nMetropolis took", ellapsed, "secs")
        plot_traces(traces, h, title="MH")
        
        # start <- Sys.time()
        # traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=1, beta=1, step=mhstep2, iter=nsamples)
        # end <- Sys.time()
        # cat("\nMetropolis rejection:",  mean(duplicated(traces[-(1:burnin),])))
        # all.traces$metropolis2 <- traces
        # all.likes$metropolis2 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        # times$metropolis2 <- end-start
        # ellapsed <- as.numeric(end-start, units="secs")
        # cat("\nMetropolis took", ellapsed, "secs")
        # 
        # start <- Sys.time()
        # traces <- sample_metropolis_h_reparam_cpp(v, W, h_init, alpha=1, beta=1, step=mhstep3, iter=nsamples)
        # end <- Sys.time()
        # cat("\nMetropolis rejection:",  mean(duplicated(traces[-(1:burnin),])))
        # all.traces$metropolis3 <- traces
        # all.likes$metropolis3 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        # times$metropolis3 <- end-start
        # ellapsed <- as.numeric(end-start, units="secs")
        # cat("\nMetropolis took", ellapsed, "secs")
        
        # MALA -----------------------------------------------------------------
        start <- Sys.time()
        traces <- sample_mala_cpp(v, W, h_init, alpha=1, beta=1, delta=0.005, iter=nsamples)
        end <- Sys.time()
        cat("\nMALA rejection:",  mean(duplicated(traces[-(1:burnin),])))
        all.traces$MALA <- traces
        all.likes$MALA <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        times$MALA <- end-start
        ellapsed <- as.numeric(end-start, units="secs")
        cat("\nMALA took", ellapsed, "secs")
        plot_traces(traces, h, title="MALA")
        
        # Gibbs ------------------------------------------------------------------------
        start <- Sys.time()
        traces <- sample_gibbs_cpp_scalable(v, W, h_init, alpha=1, beta=1, iter=nsamples)$h_n_samples
        end <- Sys.time()
        all.traces$Gibbs <- traces
        all.likes$Gibbs <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        times$Gibbs <- end-start
        ellapsed <- as.numeric(end-start, units="secs")
        cat("\nGibbs took", ellapsed, "secs")
        plot_traces(traces, h, title="Gibbs")
        
        # Hamiltonian -----------------------------------------------------------------
        start <- Sys.time()
        traces <- sample_hmc_cpp(v, W, h_init, alpha = alpha, beta = beta, epsilon = 0.0001, L=100, iter=nsamples)
        end <- Sys.time()
        cat("\nHamiltonian rejection:",  mean(duplicated(traces[-(1:burnin),])))
        all.traces$HMC <- traces
        all.likes$HMC <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        times$HMC <- end-start
        ellapsed <- as.numeric(end-start, units="secs")
        cat("\nHMC took", ellapsed, "secs")
        plot_traces(traces, h, title="HMC")

        # NUTS R ---------------------------------------------------------------
        start <- Sys.time()
        traces <- sample_nuts(v, W, theta = h_init, alpha = alpha, beta=beta, n_iter=nsamples)
        end <- Sys.time()
        cat("\n NUTS R rejection:",  mean(duplicated(traces[-(1:burnin),])))
        all.traces$nutsR <- traces
         all.likes$nutsR <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        times$nutsR <- end-start
        cat("\nNUTS R took", end-start, "secs")
        plot_traces(traces, h, title="NUTS R")
        
        # NUTS -----------------------------------------------------------------
        # start <- Sys.time()
        # traces <- sample_nuts_cpp(v, W, h_init, alpha = alpha, beta = beta, epsilon = 0.00000001, iter=nsamples)
        # end <- Sys.time()
        # cat("\n NUTS Rcpp rejection:",  mean(duplicated(traces[-(1:burnin),])))
        # all.traces$nutsCpp <- traces
        # all.likes$nutsCpp <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        # times$Cpp <- end-start
        # cat("\nHMC took", end-start, "secs")
        
        # NUTS (Stan)---------------------------------------------------------
        model <- rstan::stan_model("MMLE_H2.stan")
        #stanc("MMLE_H.stan")
        #cat(stanc("MMLE_H.stan")$cppcode)
        start <- Sys.time()
        #fit <- rstan::sampling(model, data=list(K=K, F=F, W=W, v=v, alpha=alpha, beta=beta), 
        #                       iter=nsamples, chains=1, 
        #                       verbose=TRUE, algorithm = "NUTS", init=list(list(h=h_init)))
        fit <- rstan::sampling(model, data=list(K=K, F=F, W=W, v=v, alpha=alpha, beta=beta), 
                               iter=nsamples, chains=1, 
                               verbose=TRUE, algorithm = "NUTS", init=list(list(h=h_init)),
                               control = list(stepsize = 0.01))
        end <- Sys.time()
        traces <- as.data.frame(rstan::extract(fit, permuted = FALSE, inc_warmup=TRUE)[,1,1:length(h_init)])
        traces <- as.matrix(traces)
        cat("\nNUTS rejection:",  mean(duplicated(traces[-(1:burnin),])))
        all.traces$NUTS <- traces
        all.likes$NUTS <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        times$NUTS <- end-start
        ellapsed <- as.numeric(end-start, units="secs")
        cat("\nNUTS (Stan) took", ellapsed, "secs")
        plot_traces(traces, h, title="NUTS Stan")
        
        # Hamiltonian -----------------------------------------------------------------
        # start <- Sys.time()
        # traces <- sample_hmc_cpp(v, W, h_init, alpha = alpha, beta = beta, epsilon = 0.0001, L=hmc2, iter=nsamples)
        # end <- Sys.time()
        # cat("\nHamiltonian rejection:",  mean(duplicated(traces[-(1:burnin),])))
        # all.traces$hmc2 <- traces
        # all.likes$hmc2 <- apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
        # times$hmc2 <- end-start
        # ellapsed <- as.numeric(end-start, units="secs")
        # cat("\nHMC took", ellapsed, "secs")
        # plot_traces(traces, h, title="HMC")
        
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
                           #multiESS = multiESS(mcmc(traces)),
                           RMSE    = sqrt(sum(apply(traces, 1, function(x) t(x-h)%*%(x-h)))/nsamples),
                           maxlike = max(apply(traces, 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))),
                           accepted= 1 - mean(duplicated(traces))) 
          
          df.summary.last <- rbind(df.summary.last, df)
        }
        
        
        # Store all traces ready to be plotted
        df.traces <- data.frame()
        for(i in 1:length(all.traces)){
          method <- names(all.traces)[i]
          traces <- all.traces[[i]]
          nsamples_ <- nrow(traces)
          milsecs <- seq(0,as.double(times[[method]]), length.out=nsamples_)*1000
          df.traces <- rbind(df.traces, 
                             as.data.frame(traces) %>%
                             mutate(method     = method, 
                                    iteration  = 1:nsamples_, 
                                    time=milsecs) %>%
                               gather(k, value, -method, -iteration, -time))
        }
        
        # Plot traces of all methods
        p <- ggplot(df.traces, aes(x=iteration, y=value, color=k)) + geom_line() + 
            geom_hline(yintercept = h, linetype='dashed') +
            facet_grid(method~.) +
            #scale_y_log10()+
            theme_bw() + theme(legend.position = "none")
        print(p)
        
        
        # Store Likelihoods (not used) -------------------------------------------------
        res <- list()
        for(method in names(all.traces)){
          likelihoods <- apply(all.traces[[method]], 1, function(x) posterior_h_cpp_(x, v, W, alpha, beta))
          res[[method]] <- data.frame(method = method, sample=1:length(likelihoods), likelihood = likelihoods )
        }
        

        
      }
    }
  }
} ## end loop

# Metric results for paper
if(FALSE){
df.summary <- df.summary.last
df.summary$method <- as.character(df.summary$method)
cat("\n\nEnd of experiments. Plotting the results")
df.summary.backup <- df.summary #just in case
df.summary$time <- as.numeric(df.summary$time)

df.summary = df.summary %>% mutate(K= paste0("K = ", K))

#df.summary[df.summary$method == "metropolis1",]$method <- "MH-1"
#df.summary[df.summary$method == "metropolis2",]$method <- "MH-2"
#df.summary[df.summary$method == "metropolis3",]$method <- "MH-3"
#df.summary[df.summary$method == "gibbs",]$method <- "Gibbs"
#df.summary[df.summary$method == "hmc1",]$method <- "HMC-1"
#df.summary[df.summary$method == "hmc2",]$method <- "HMC-2"
#df.summary[df.summary$method == "nuts",]$method <- "NUTS"
#df.summary[df.summary$method == "vb",]$method <- "VB"


df.summary$K <- as.factor(df.summary$K)
df.summary$F <- as.factor(df.summary$F)
df.summary$method <- as.factor(df.summary$method)

#levels <- c("VB", "MH-1","MH-2","MH-3","HMC-1","HMC-2","Gibbs","NUTS")
#df.summary$method <- factor(df.summary$method, levels = levels)

mytheme <- theme_bw() + theme(aspect.ratio=1/2.5, 
                              axis.title.y=element_blank(),
                              legend.title = element_blank(), 
                              legend.position = "none",
                              strip.background = element_blank(),
                              text = element_text(size=20),
                              plot.title = element_text(hjust = 0.5))
                # plot.margin=grid::unit(c(0,0,0,0), "mm"))

mytheme_legend <- theme_bw() + theme(aspect.ratio=1/2.5, 
                                     axis.title.y=element_blank(),
                                     legend.title = element_blank(),
                                     strip.background = element_blank(),
                                     text = element_text(size=20), 
                                     plot.title = element_text(hjust = 0.5),
                                     legend.position = "bottom")

mywidth=15
myheight=15
# RMSE #######
# w.r.t F
p1 <- ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=RMSE, fill=method, color=method)) +  
  facet_grid(K~.) +
  #scale_y_log10() +
  scale_shape_manual(values=1:nlevels(df.summary$method)) +
  #scale_color_manual(values=1:nlevels(df.summary$method), guide=FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) + 
  mytheme_legend + ggtitle("RMSE") 

p1                     
ggsave("RMSE_synthetic_hard_col.eps",
       scale = 1, units = "cm", dpi = 300, width=mywidth, height=myheight)

# worst_z #######
p2 <- ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=worst_z, fill=method, color=method)) +  
  facet_grid(K~.) +
  scale_y_log10() +
  scale_shape_manual(values=1:nlevels(df.summary$method)) +
  #scale_color_manual(values=1:nlevels(df.summary$method), guide=FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) + 
  #ylim(0,50) +
  geom_hline(yintercept = c(1.96), linetype='dashed') +
  mytheme_legend + ggtitle("max. Geweke z-score")
p2
ggsave("synthetic_hard_worst_z_col.eps",
       scale = 1, units = "cm", dpi = 300, width=mywidth, height=myheight)


# ESS/time ########
p3 <- ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=minESS/time, fill=method, color=method)) +  
  facet_grid(K~.) +
  scale_shape_manual(values=1:nlevels(df.summary$method)) +
  #scale_color_manual(values=1:nlevels(df.summary$method), guide=FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) +
  scale_y_log10() +
  mytheme_legend + ggtitle("min. ESS/time")


ggsave("synthetic_hard_esstime_col.eps",
       scale = 1, units = "cm", dpi = 300, width=mywidth, height=myheight)

# ESS ############
# w.r.t F
p4 <- ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=minESS, fill=method, color=method)) +  
  facet_grid(K~.) +
  scale_y_log10()+
  scale_shape_manual(values=1:nlevels(df.summary$method)) +
  #scale_color_manual(values=1:nlevels(df.summary$method), guide=FALSE) +
  scale_colour_brewer(palette = "Dark2") +
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) + 
  mytheme_legend + ggtitle("min. ESS")

ggsave("synthetic_hard_ess_col.eps",
       scale = 1, units = "cm", dpi = 300, width=mywidth, height=myheight)

pall <- grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 4, nrow = 1)
ggplot2::ggsave("synthetic_hard_all.eps", pall,
       scale = 1, units = "cm", dpi = 300, width=50, height=20) #50x20
###################
#*******************************************************************************
# NOT USED IN THE PAPER
#*******************************************************************************


  
# maxlikelihood #######
ggplot(df.summary %>% filter(alphadir==alphadir), 
       aes(x=F, y=maxlike, fill=method)) +  
  facet_grid(K~.) +
  scale_y_log10()+
  geom_point(position = position_dodge(width=0.75), aes(group=method, shape=method)) + 
  mytheme
}
