//
// Stan model for logistic regression
//
data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  int<lower=0,upper=1> y[N];
  vector[p] mu_beta;
  cov_matrix[p] Sigma_beta;
}

parameters {
  vector[p] beta;
}

model {
  y ~ bernoulli_logit(X*beta);
  beta ~ multi_normal(mu_beta, Sigma_beta);
}
