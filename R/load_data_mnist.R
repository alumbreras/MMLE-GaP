# Load the MNIST digit recognition dataset into R
# http://yann.lecun.com/exdb/mnist/
# assume you have all 4 files and gunzip'd them
# creates train$n, train$x, train$y  and test$n, test$x, test$y
# e.g. train$x is a 60000 x 784 matrix, each row is one digit (28x28)
# call:  show_digit(train$x[5,])   to see a digit.
# brendan o'connor - gist.github.com/39760 - anyall.org

show_digit_old <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}

show_digit <- function(arr784, col=gray(255:1/255), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}


load_mnist <-function(){
  filename <- '../data/mnist/train-images.idx3-ubyte'
  f = file(filename,'rb')
  readBin(f,'integer', n=1, size=4, endian='big')
  n = readBin(f,'integer', n=1, size=4, endian='big')
  nrow = readBin(f, 'integer', n=1, size=4, endian='big')
  ncol = readBin(f, 'integer', n=1, size=4, endian='big')
  x = readBin(f, 'integer', n=n*nrow*ncol, size=1, signed=FALSE)
  x = matrix(x, ncol=nrow*ncol, byrow=TRUE)
  close(f)
  return(x)
}


if(FALSE){
  par(pty="s", mfrow=c(2,2))
  show_digit(x[1,])
  show_digit(x[2,])
  show_digit(x[3,])
  show_digit(x[40,])
}