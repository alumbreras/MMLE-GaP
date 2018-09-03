
#' Loads the data matrix (V matrix) from the training file
#' @return A matrix with one row per image
#' @export
load_image_file <- function(filename){
  f = file(filename, "rb")
  magicn = readBin(f, 'integer', n=1, size=4, endian='big') # magic number
  nitems = readBin(f, 'integer', n=1, size=4, endian='big') # number of items
  nrows = readBin(f,'integer',n=1,size=4,endian='big')
  ncols = readBin(f,'integer',n=1,size=4,endian='big')
  V <- do.call('cbind', lapply(1:nitems, function(i) { 
                        readBin(f,'integer',n=nrows*ncols,size=1,signed=F)
                        }))
  close(f)
  return(V)
}


#' Plot an image from a given row in the data matrix
#' @param V the whole data matrix
#' @param n the number of the column to plot
#' @details It gets the n-th column, loads it to a 28x28 matrix and then plots it
plot_image <- function(data, n){
  x <- data[,n]
  m <- matrix(x, 28, 28, byrow = TRUE)
  image(t(m)[,28:1], col=gray(255:1/255))
}

# Example of PCA reconstruction with eigennumbers
# see:
# http://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
if(FALSE){

  train_images_file <- "train-images.idx3-ubyte"
  train_labels_file <- "train-labels.idx1-ubyte"
  test_images_file <- "t10k-images.idx3-ubyte"
  test_labels_file <- "t10k-labels.idx1-ubyte"
  
  V <- load_image_file(train_images_file)
  par(mfrow=c(3,3))
  for(i in 1:9){
   plot_image(V,i)
  }
  
  # Data input. Data points in columns, rows are dimensions (or features)
  V <- V[,1:10000]
  
  # Compute PCA
  pca <- prcomp(t(V)) # the function needs features columns
  
  # Rebuild the original data using only some of the dimensions
  ncomp <- 250
  W <- pca$rotation[,1:ncomp] # dictionary, eigenvectors U in PCA notation
  Y <- t(pca$x[,1:ncomp]) # projections, Y = U^T X
  mu = rowMeans(V) # data mean to reconstruct X properly
  Xhat = W %*% Y 
  Xhat = t(scale(t(Xhat), center = -mu, scale = FALSE))
  
  par(mfrow=c(3,3))
  for(i in 1:9){
    plot_image(W,i)
  }
  
  # Plot reconstructions
  par(mfrow=c(2,2))
  for(i in 1:10){
    plot_image(V,i)
    plot_image(Xhat,i)
  }
}


if(FALSE){
  data <- load_synthetic()
  save(data, file = 'synthetic_V_W_F.RData')
}