# ************************************************************
# Example of MMLE - NMF exploiting the closeform expression
# Plots 2-D traces of the likelihood of W
# ************************************************************
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(fields)
library(rgl)

set.seed(1)

################################## DATA ########################################
# Generate syhtetic data ####
data <- load_synthetic()
V <- data$V
#V <- V[1:4, 1:10]
F <- dim(V)[1]
N <- dim(V)[2]
plot_multiple_W(data$V, data$W, data$H, cols=3)


# generate synthetic data from a GaP model ####
F <- 1
K <- 2
V <- matrix(c(1,1,10,10), nrow=F)
initW <- matrix(1, nrow = F, ncol = K) # K initially overestimated

alpha <- 1

maxiters = 1000
nsamples = 250
################################## ALGOS #######################################
maxiters = 100
nsamples = 12
# METROPOLIS ####
res.mh <- nmf_mmle_mcem(V, initW, maxiters, nsamples, 
                        algo='metropolis', 
                        mincrease = "constant", 
                        alpha=alpha)
res <- res.mh
samples <- res$traces$W

# Plot result matrix
plot_multiple_W(W, initW, res$W)
plot(sort(colSums(res$W), decreasing = TRUE)) # plot column norms

# Plot traces
df <- as.data.frame(samples)
df$sample <- 1:nrow(df)
df <- gather(df, k, norm, -sample)
df <- df[complete.cases(df),]
ggplot(df, aes(x=sample, y=norm, color=k)) + geom_line()

# Plot samples over the true function
Ngrid <-100
w11 <- seq(0,10,length.out=Ngrid)
w12 <- seq(0,10,length.out=Ngrid)
xygrid <- expand.grid(x=w11, y=w12)
df.xyz <- xygrid %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(x,y), nrow=1, ncol=2), 
                            alpha=alpha, beta=1, log= FALSE, vectorizable=TRUE)) %>%
  ungroup

# Plot true likelihood p(V | W)
z <- df.xyz$z
z <- matrix(z, nrow=sqrt(length(z)))
image.plot(z)
contour(z, add = TRUE)

# Plot samples of W
par(pty="s") # squared plots
contour(w11, w12, z, xlim=c(0,100), ylim=c(0,100))
contour(w11, w12, z)
lines(samples, col='yellow')
points(samples, col='red', pch=19, cex=0.2)

# GIBBS ####
res.gibbs <- nmf_mmle_mcem(V, initW, maxiters, nsamples, 
                           algo='gibbs', 
                           mincrease = "constant", 
                           alpha=alpha)
res <- res.gibbs
samples <- res$traces$W

# Plot result matrix
plot_multiple_W(W, initW, res$W)
plot(sort(colSums(res$W), decreasing = TRUE)) # plot column norms

# Plot traces
df <- as.data.frame(samples)
df$sample <- 1:nrow(df)
df <- gather(df, k, norm, -sample)
df <- df[complete.cases(df),]
ggplot(df, aes(x=sample, y=norm, color=k)) + geom_line()

# Plot samples over the true function
Ngrid <-100
w11 <- seq(0,10,length.out=Ngrid)
w12 <- seq(0,10,length.out=Ngrid)
xygrid <- expand.grid(x=w11, y=w12)
df.xyz <- xygrid %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(x,y), nrow=1, ncol=2), 
                            alpha=alpha, beta=1, log= FALSE, vectorizable=TRUE)) %>%
  ungroup

z <- df.xyz$z
z <- matrix(z, nrow=sqrt(length(z)))
image.plot(z)
contour(z, add = TRUE)

par(pty="s") # squared plots
contour(w11, w12, z, xlim=c(0,100), ylim=c(0,100))
contour(w11, w12, z)
lines(samples, col='yellow')
points(samples, col='red', pch=19, cex=0.2)



#############################""
# HAMILTONIAN ####
steps <- array(NA, dim=c(N*L, K))
nn <- 1
proposals <- array(NA, dim=c(N, K))
jj <- 1

# For F=1 peta en stan por alguna chorrada de dimensiones.
# Usar uno de los datos sinteticos de F=2 o mas
res.hmc <- nmf_mmle_mcem(V, initW, maxiters=100, nsamples=550, 
                         algo='hamiltonian-stan', 
                         mincrease = "constant", 
                         alpha=alpha,
                         L=5, epsilon=0.5)
res <- res.hmc
samples <- res$traces$W

# Plot result matrix
plot_multiple_W(W, initW, res$W)
plot(sort(colSums(res$W), decreasing = TRUE)) # plot column norms

# Plot traces
df <- as.data.frame(samples)
df$sample <- 1:nrow(df)
df <- gather(df, k, norm, -sample)
df <- df[complete.cases(df),]
ggplot(df, aes(x=sample, y=norm, color=k)) + geom_line()

# Plot samples over the true function
Ngrid <-100
w11 <- seq(0,10,length.out=Ngrid)
w12 <- seq(0,10,length.out=Ngrid)
xygrid <- expand.grid(x=w11, y=w12)
df.xyz <- xygrid %>% 
  rowwise %>% 
  mutate(z = likelihood_V_W(V,  W = matrix(c(x,y), nrow=1, ncol=2), 
                            alpha=alpha, beta=1, log= FALSE, vectorizable=TRUE)) %>%
  ungroup

# Plot true likelihood p(V | W)
z <- df.xyz$z
z <- matrix(z, nrow=sqrt(length(z)))
image.plot(z)
contour(z, add = TRUE)

nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(z, nbcol)
persp3d(w11, w12, z, col=color[zcol])

# Plot samples of W
par(pty="s") # squared plots
contour(w11, w12, z, xlim=c(0,100), ylim=c(0,100))
contour(w11, w12, z)
lines(samples, col='yellow')
points(samples, col='red', pch=19, cex=0.2)