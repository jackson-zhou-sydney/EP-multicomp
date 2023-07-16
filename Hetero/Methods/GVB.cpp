#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nljl(vec theta, mat X_1, mat X_2, vec y, mat Sigma_theta, vec mu_theta) {
  // Negative log joint likelihood
  int p_1 = X_1.n_cols;
  int p_2 = X_2.n_cols;
  
  vec beta_1 = theta.head(p_1);
  vec beta_2 = theta.tail(p_2);
  
  vec y_XB_1 = y - X_1*beta_1;
  vec XB_2 = 2.0*X_2*beta_2;
  vec theta_mu_theta = theta - mu_theta;
  
  return 0.5*(accu(XB_2) + dot(y_XB_1, y_XB_1/exp(XB_2)) + dot(theta_mu_theta, inv_sympd(Sigma_theta)*theta_mu_theta));
}

// [[Rcpp::export]]
vec nljl_grad(vec theta, mat X_1, mat X_2, vec y, mat Sigma_theta, vec mu_theta) {
  // Gradient of negative log joint likelihood
  int n = X_1.n_rows;
  int p_1 = X_1.n_cols;
  int p_2 = X_2.n_cols;
  
  vec beta_1 = theta.head(p_1);
  vec beta_2 = theta.tail(p_2);
  
  vec y_XB_1 = y - X_1*beta_1;
  vec XB_2 = 2.0*X_2*beta_2;
  vec y_XB_XB = y_XB_1/exp(XB_2);
  
  vec a_1 = X_1.t()*y_XB_XB;
  mat X_3 = zeros(n, p_2);
  for (int j = 0; j < p_2; ++j) {
    X_3.col(j) = 2.0*X_2.col(j)%y_XB_XB;
  }
  vec a_2 = 0.5*(X_3.t()*y_XB_1 - sum(2.0*X_2.t(), 1));
  
  return inv_sympd(Sigma_theta)*theta - join_cols(a_1, a_2);
}
