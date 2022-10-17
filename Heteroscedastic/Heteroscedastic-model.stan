//
// Stan model for heteroscedastic linear regression
//
data {
  int<lower=0> N;
  int<lower=0> p_1;
  int<lower=0> p_2;
  matrix[N, p_1] X_1;
  matrix[N, p_2] X_2;
  vector[N] y;
  vector[p_1 + p_2] mu_theta;
  cov_matrix[p_1 + p_2] Sigma_theta;
}

parameters {
  vector[p_1 + p_2] theta;
}

transformed parameters {
  vector[p_1] beta_1;
  vector[p_2] beta_2;
  beta_1 = theta[1:p_1];
  beta_2 = theta[(p_1 + 1):(p_1 + p_2)];
}

model {
  y ~ normal(X_1*beta_1, exp(X_2*beta_2));
  theta ~ multi_normal(mu_theta, Sigma_theta);
}
