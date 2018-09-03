# Compares the R and Rcpp implementations of the likelihood

V <- matrix(c(0,0,0,1,1,
              0,4,8,10,50,
              0,6,9,30,100), nrow=3, byrow = TRUE)

C <- matrix(c(1,10,
              1,0,
              1,0), nrow=3, byrow = TRUE)

W <- matrix(c(0,0,10,
              0,4,10,
              0,6,100), nrow=3, byrow = TRUE)


#likelihood_C_W(C, W, alpha=1, beta=1, log=TRUE)
like = likelihood_V_W(V, W, alpha=2, beta=3, log=TRUE)
cat("\nlikelihood V W: ", like)

like = likelihood_V_W_cpp(V, W,  alpha=2, beta=3)
cat("\nlikelihood V W C++: ", like)

rbenchmark::benchmark(likelihood_V_W(V, W, alpha=1, beta=1, log=TRUE),  
                      likelihood_V_W_cpp(V, W), alpha=1, beta=1)