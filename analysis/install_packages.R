# Install packages that may be used in the scripts under the /analytics folder
list.of.packages <- c("dplyr", "tidyr", "reshape2", "ggplot2", "cowplot",
						"devtools", "profvis", "Rcpp", "RcppArmadillo", "RcppEigen")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


