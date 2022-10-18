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
  vector[p + 2] mu_bk;
  cov_matrix[p + 2] Sigma_bk;
}

parameters {
  vector[p + q + 2] theta;
}

transformed parameters {
  vector[p] beta;
  vector[q] u;
  real kappa_y;
  real kappa_u;
  vector[p + 2] bk;
  
  beta = theta[1:p];
  u = theta[(p + 1):(p + q)];
  kappa_y = theta[p + q + 1];
  kappa_u = theta[p + q + 2];
  bk = append_row(append_row(beta, kappa_y), kappa_u);
}

model {
  y ~ normal(X*beta + Z*u, rep_vector(exp(kappa_y), N));
  u ~ normal(rep_vector(0, q), rep_vector(exp(kappa_u), q));
  bk ~ multi_normal(mu_bk, Sigma_bk);
}
