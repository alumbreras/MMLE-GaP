data {
  int K;
  int F;
  matrix[F,K] W;
  int v[F];
  real alpha;
  real beta;
}

transformed data {
	vector[K] alphas;
	vector[K] betas;
	for(k in 1:K){
		alphas[k] = alpha;
		betas[k] = beta;
	}
}

parameters {
  vector<lower=0>[K] h;
}

transformed parameters {
	vector[F] lambdas;
	lambdas = W * h;
}

model {    
    h ~ gamma(alphas, betas); // shape, rate
    v ~ poisson(lambdas);
}
