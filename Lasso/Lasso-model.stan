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

model {
  y ~ normal_id_glm(X, 0, theta[1:p], exp(theta[p + 1]));
  theta[1:p] ~ double_exponential(0, exp(theta[p + 1])/lambda);
  theta[p + 1] ~ normal(mu_kappa, sqrt(sigma_2_kappa));
}
