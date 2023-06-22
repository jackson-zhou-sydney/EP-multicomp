//
// Stan model for quantile regression
//
functions {
  vector rho(vector y, real tau) {
    return 0.5*(fabs(y) + (2*tau - 1)*y);
  }
}

data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N, p] X;
  vector[N] y;
  vector[p + 1] mu_theta;
  cov_matrix[p + 1] Sigma_theta;
  real<lower=0,upper=1> tau;
}

parameters {
  vector[p + 1] theta;
}

model {
  target += -N*theta[p + 1] - sum(rho(y - X*theta[1:p], tau))/exp(theta[p + 1]);
  theta ~ multi_normal(mu_theta, Sigma_theta);
}
