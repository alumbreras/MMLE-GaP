show_digit <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}

url <- "http://yann.lecun.com/exdb/mnist/train-images-idx3-ubyte.gz"
filename <- basename(url)
download.file(url = url, destfile = filename)
f <- gzfile(filename, "rb")
readBin(f, integer(), n = 1, size = 4, endian = "big")
n <- readBin(f, integer(), n = 1, size = 4, endian = "big")
nrow <- readBin(f, integer(), n = 1, size = 4, endian = "big")
ncol <- readBin(f, integer(), n = 1, size = 4, endian = "big")
x <- readBin(f, integer(), n = n * nrow * ncol, size = 1, signed = FALSE)
x <- matrix(x, ncol = nrow * ncol, byrow = TRUE)
close(f)

show_digit(x[1,])


