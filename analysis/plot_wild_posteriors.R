# Draws random h, W, v and plots the posterior p(h | ) 
# The goal is to find all kind of weird posteriors

library(gtools)
library(gtools)
library(dplyr)
library(tidyr)
library(fields)
library(rstan)
library(coda)
library(fields)
library(rgl)
library(gtools)
library(colorRamps)

alpha <- 1
beta <- 1
K <- 2
F <- 10
alphadir <- 100 # high: uniform vectors # low (<1): sparse vectors
#alphadir <- 1 # high: uniform vectors # low (<1): sparse vectors


datum <- list()

for(i in 1:10){
# Data 1-------------------------------------------------------------------------


#case 4
W <- abs(matrix(rnorm(F*K, 10, 100), ncol=K))

# case 3
W <- matrix(runif(F*K), ncol=K)

# case 3
W <- matrix(rexp(F*K, 0.1),  ncol=K)

# synthetic 
W <- matrix(c(1,0,0,1), nrow=2)
W <- matrix(c(1,1,0,0), nrow=2)
W <- matrix(c(1,0,1,0), nrow=2)


# case 3
W <- matrix(runif(F*K), ncol=K)

# case 1
W <- rdirichlet(F, alpha = rep(alphadir,K))
W <- W*1

# case 2
W <- t(rdirichlet(K, alpha = rep(alphadir,F)))

# Generate observation v
W <- W*1
h <- rgamma(K, shape=100, rate=1)
v <- rpois(F, lambda = W%*%h)

datum[[i]] <- list(W=W, v=v, h=h)
}

par(pty="s", mfrow=c(3,3))
plot.new()
for(i in 1:length(datum)){
  W <- datum[[i]]$W
  v <- datum[[i]]$v
  
  Ngrid <-100
  valmax <- 200
  h1 <- seq(0,valmax,length.out=Ngrid)
  h2 <- seq(0,valmax,length.out=Ngrid)
  xygrid <- expand.grid(x=h1, y=h2)
  df.xyz <- xygrid %>% 
    rowwise %>% 
    mutate(z = posterior_h_cpp_(h = c(x,y), v=v, W=W, alpha=alpha, beta=1)) %>%
    #mutate(z = posterior_h(h = c(x,y), v=v, W=W, alpha=alpha, beta=1, log=TRUE)) %>%
    ungroup
  z <- df.xyz$z-max(df.xyz$z)
  z <- exp(df.xyz$z-max(df.xyz$z))
  z <- matrix(z, nrow=sqrt(length(z)))
  
  # Plot on screen
  #par(pty="s")
  #plot.new()
  #image.plot(z)
  #contour(h1, h2, z, add = TRUE)

  
  # Plot on screen
  #filled.contour(h1, h2, z, xlab="h[1]", ylab="h[2]", key.axes = NA, color.palette = matlab.like)
  contour(h1, h2, z, add = FALSE)
  
}





if(FALSE){
  z <- df.xyz$z-max(df.xyz$z)
  z <- exp(df.xyz$z-max(df.xyz$z))
  z <- matrix(z, nrow=sqrt(length(z)))
  nbcol = 100
  color = rev(rainbow(nbcol, start = 0/6, end = 4/6))# Samplers ---------------------------------------------------------------------
  zcol  = cut(z, nbcol)
  persp3d(h1, h2, z, col=color[zcol])
}
