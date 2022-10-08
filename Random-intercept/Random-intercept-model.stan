//
// Stan model for random intercept linear regression
//
data {
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> q;
  matrix[N, p] X;
  matrix[N, q] Z;
  vector[N] y;
  vector[p] mu_beta;
  real mu_tau_y;
  real mu_tau_u;
  cov_matrix[p] Sigma_beta;
  real<lower=0> sigma_2_tau_y;
  real<lower=0> sigma_2_tau_u;
}

parameters {
  vector[p + q + 2] theta;
}

transformed parameters {
  vector[p] beta;
  vector[q] u;
  real tau_y;
  real tau_u;
  
  beta = theta[1:p];
  u = theta[(p + 1):(p + q)];
  tau_y = theta[p + q + 1];
  tau_u = theta[p + q + 2];
}

model {
  y ~ normal(X*beta + Z*u, rep_vector(exp(tau_y), N));
  beta ~ multi_normal(mu_beta, Sigma_beta);
  u ~ normal(rep_vector(0, q), rep_vector(exp(tau_u), q));
  tau_y ~ normal(mu_tau_y, sqrt(sigma_2_tau_y));
  tau_u ~ normal(mu_tau_u, sqrt(sigma_2_tau_u));
}
