data {
  int K;
  int F;
  matrix[F,K] W;
  int v[F];
}
parameters {
  real<lower=0> h[K];
}
model {
	real lambda;    
    for (k in 1:K) {
		h[k] ~ gamma(1,1);
    }
    
	for (f in 1:F) {
		lambda = 0;
		for (k in 1:K) {
			lambda = lambda + W[f,k] * h[k]; 
		}
		v[f] ~ poisson(lambda);
		}
}
