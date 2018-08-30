library(rgl)
# Generate N samples of v (F=1, K=1). 
# Since v (F=1) is the sum of K Negative Binomials, for K=1 it is just a 
# Negative Binomial
N  <- 100
Nx <- 100
alpha_star <- 10
alpha <- 50
beta       <- 1
w_star     <- 0.2

v <- rnbinom(N, size=alpha_star, prob=beta/(w_star+beta))
V <- t(as.matrix(v))
hist(v, breaks=100)

# Only through a line
if(TRUE){
M <- mean(v)/alpha
w12 <- seq(0,M/2,length.out=Nx)
w11 <- M-w12

xygrid <- data.frame(w12, w11) %>% filter(w11>=0, w12>=0)
df.xyz <- xygrid %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(x,y), nrow=1, ncol=2), 
                            alpha=alpha, beta=beta, log=TRUE, vectorizable=TRUE)) %>%
  ungroup

plot(df.xyz$w11, df.xyz$z)
}

if(FALSE){
# Compute the MLE of w by grid search, assuming K=2
# The goal is to see whether the MLE of the estimator w=(w1,w2) has one 
# of its components at zero.
w11 <- seq(0,1,length.out=Nx)
w12 <- seq(0,1,length.out=Nx)
xygrid <- expand.grid(x=w11, y=w12)

df.xyz <- xygrid %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(x,y), nrow=1, ncol=2), 
                            alpha=alpha, beta=beta, log=FALSE, vectorizable=TRUE)) %>%
  ungroup



df.xyz %>% arrange(-z)
MLE_probs <- df.xyz[which.max(df.xyz$z),]

z <- df.xyz$z
z <- matrix(z, nrow=sqrt(length(z)))
image.plot(z)
contour(z, add = TRUE)
title(MLE_probs[1:2])

zcol  = cut(z, nbcol)
persp3d(w11, w12, z, col=color[zcol])
}