//
// Stan model for probit regression
//
data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  int<lower=0,upper=1> y[N];
  vector<lower=0,upper=0>[p] mu_beta;
  cov_matrix[p] Sigma_beta;
}

parameters {
  vector[p] beta;
}

model {
  y ~ bernoulli(Phi(X*beta));
  theta ~ multi_normal(mu_beta, Sigma_beta);
}
