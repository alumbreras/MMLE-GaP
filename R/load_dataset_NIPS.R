# Load the NIPS dataset
load_nips <-function(){
  library(reshape2)
  library(Matrix)
  
  # Read the vocabulary. This is useless for now.
  # It will be useful when analyzing the results
  filename <- '../data/nips/vocab.nips.txt'
  vocabulary <- readLines(filename)

  filename <- '../data/nips/docword.nips_noheader.txt'
  dat <- read.table(filename)
  names(dat) <- c("document", "word", "count")
  V = sparseMatrix(dat$word, dat$document, x = dat$count)
  # V <- acast(dat, word~document, value.var='count', fill=0) # Not practical for large data

  return(V)
}


