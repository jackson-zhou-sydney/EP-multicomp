//
// Stan model for quantile regression
//
functions {
  real rho(real y, real tau) {
    return 0.5*(fabs(y) + (2*tau - 1)*y);
  }
  real asym_double_exponential_lpdf(real y, real mu, real sigma, real tau) {
    return log(tau) + log(1 - tau) - log(sigma) - rho((y - mu)/sigma, tau);
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

transformed parameters {
  vector[p] beta;
  real kappa;
  
  beta = theta[1:p];
  kappa = theta[p + 1];
}

model {
  for (i in 1:N) {
    y[i] ~ asym_double_exponential(X[i]*beta, exp(kappa), tau);
  }
  theta ~ multi_normal(mu_theta, Sigma_theta);
}
