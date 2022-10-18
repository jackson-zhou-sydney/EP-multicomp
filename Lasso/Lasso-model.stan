//
// Stan model for lasso linear regression
//
data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  vector[N] y;
  real mu_kappa;
  real<lower=0> sigma_2_kappa;
  real<lower=0> lambda;
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
  beta ~ double_exponential(rep_vector(0, p), rep_vector(exp(kappa)/lambda, p));
  kappa ~ normal(mu_kappa, sqrt(sigma_2_kappa));
}
