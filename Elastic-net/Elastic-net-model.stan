//
// Stan model for elastic net linear regression
//
data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  vector[N] y;
  real mu_kappa;
  real<lower=0> sigma_2_kappa;
  real<lower=0> lambda_1;
  real<lower=0> lambda_2;
}

parameters {
  vector[p + 1] theta;
}

transformed parameters {
  vector[p] beta;
  real kappa;
  
  beta = theta[1:p];
  kappa = theta[p + 1];
}

model {
  y ~ normal(X*beta, rep_vector(exp(kappa), N));
  beta ~ double_exponential(rep_vector(0, p), rep_vector(exp(kappa)/lambda_1, p));
  beta ~ normal(rep_vector(0, p), rep_vector(sqrt(0.5*exp(kappa)/lambda_2), p));
  kappa ~ normal(mu_kappa, sqrt(sigma_2_kappa));
}
