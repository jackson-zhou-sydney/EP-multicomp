//
// Stan model for elastic net linear regression
//
data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  vector[N] y;
  real kappa_1;
  real kappa_2;
  real mu_tau;
  real<lower=0> sigma_2_tau;
}

parameters {
  vector[p + 1] theta;
}

transformed parameters {
  vector[p] beta;
  real tau;
  
  beta = theta[1:p];
  tau = theta[p + 1];
}

model {
  y ~ normal(X*beta, rep_vector(exp(tau), N));
  beta ~ double_exponential(rep_vector(0, p), rep_vector(exp(tau - kappa_1), p));
  beta ~ normal(rep_vector(0, p), rep_vector(sqrt(0.5*exp(tau - kappa_2)), p));
  tau ~ normal(mu_tau, sqrt(sigma_2_tau));
}
