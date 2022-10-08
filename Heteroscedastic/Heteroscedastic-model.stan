//
// Stan model for heteroscedastic linear regression
//
data {
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> q;
  matrix[N, p] X;
  matrix[N, q] Z;
  vector[N] y;
  vector[p + q] mu_theta;
  cov_matrix[p + q] Sigma_theta;
}

parameters {
  vector[p + q] theta;
}

transformed parameters {
  vector[N] latent_beta;
  vector[N] latent_alpha;
  latent_beta = X*theta[1:p];
  latent_alpha = Z*theta[(p + 1):(p + q)];
}

model {
  y ~ normal(latent_beta, sqrt(exp(latent_alpha)));
  theta ~ multi_normal(mu_theta, Sigma_theta);
}
