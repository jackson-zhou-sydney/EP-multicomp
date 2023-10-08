#include <RcppArmadillo.h>
#include <RcppNumerical.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace arma;
using namespace Eigen;
using namespace Numer;
using namespace std;

typedef Map<MatrixXd> MapMat;
typedef Map<VectorXd> MapVec;

class Neg_Expnt: public MFuncGrad {
private:
  double A;
  double B;
  double C;
  double D;
public:
  Neg_Expnt(double A_, double B_, double C_, double D_) : A(A_), B(B_), C(C_), D(D_) {}
  
  double f_grad(Constvec& x, Refvec grad) {
    double x_dbl = x(0);
    double f = A*x_dbl + pow(x_dbl - B, 2.0)/(2.0*C) + D/(2.0*exp(2.0*x_dbl));
    grad(0) = A + (x_dbl - B)/C - D/exp(2.0*x_dbl);
    return f;
  }
};

double trap_unif(double delta_x, rowvec y) {
  // Trapezoidal integration with uniform intervals
  y(0) = 0.5*y(0);
  y.back() = 0.5*y.back();
  return delta_x*sum(y);
}

double E_lnig(double A, double B, double C, double D, String fun, int n_grid) {
  // Expectation of log-normal and inverse gamma product
  Neg_Expnt neg_expnt(A, B, C, D);
  
  VectorXd x_est = VectorXd::Zero(1);
  double fopt;
  int status = optim_lbfgs(neg_expnt, x_est, fopt, 1000);
  if (status < 0) {
    stop("Failed to converge");
  }
  
  double mu = x_est(0);
  double sigma = 1.0/sqrt((1.0/C + 2*D*exp(-2.0*mu)));
  
  double lb = mu - 5.0*sigma;
  double ub = mu + 5.0*sigma;
  
  vec x_values = linspace(lb, ub, n_grid);
  double delta_x = x_values(1) - x_values(0);
  mat y_matrix = zeros(2, n_grid);
  
  if (fun == "x") {
    for (int i = 0; i < n_grid; ++i) {
      double x = x_values(i);
      y_matrix(0, i) = exp(-A*x - pow(x - B, 2.0)/(2.0*C) - D/(2.0*exp(2.0*x)) + fopt);
      y_matrix(1, i) = x*y_matrix(0, i);
    }
  } else if (fun == "x^2") {
    for (int i = 0; i < n_grid; ++i) {
      double x = x_values(i);
      y_matrix(0, i) = exp(-A*x - pow(x - B, 2.0)/(2.0*C) - D/(2.0*exp(2.0*x)) + fopt);
      y_matrix(1, i) = pow(x, 2.0)*y_matrix(0, i);
    }
  } else if (fun == "1/exp(2*x)") {
    for (int i = 0; i < n_grid; ++i) {
      double x = x_values(i);
      y_matrix(0, i) = exp(-A*x - pow(x - B, 2.0)/(2.0*C) - D/(2.0*exp(2.0*x)) + fopt);
      y_matrix(1, i) = exp(-2.0*x)*y_matrix(0, i);
    }
  } else {
    stop("fun must be one of: x, x^2, or 1/exp(2*x)");
  }
  
  return trap_unif(delta_x, y_matrix.row(1))/trap_unif(delta_x, y_matrix.row(0));
}

// [[Rcpp::export]]
List mfvb(mat X, vec y, double sigma_2_kappa, double mu_kappa,
          double lambda, int min_iter, int max_iter, double thresh, int n_grid, bool verbose) {
  // MFVB for Bayesian lasso linear regression
  int n = X.n_rows;
  int p = X.n_cols;
  
  double bd_mu;
  double bd_Sigma;
  double bd_a;
  
  mat XTX = X.t()*X;
  vec XTy = X.t()*y;
  
  // Initialisations
  vec mu_beta = zeros(p);
  mat Sigma_beta = eye(p, p);
  double E_ie_2_kappa = 1;
  vec E_a = ones(p);
  mat Q;
  mat Q_inv;
  
  bool converged = false;
  
  // Main MFVB loop
  for (int i = 0; i < max_iter; ++i) {
    // Store old values
    vec mu_beta_old = mu_beta;
    mat Sigma_beta_old = Sigma_beta;
    vec E_a_old = E_a;
    
    if (verbose) {
      Rcout << "---- Current iteration: " << i << " ----\n";
    }
    
    // Update q(beta)
    Q_inv = XTX + pow(lambda, 2.0)*diagmat(E_a);
    Q = inv(Q_inv);
    mu_beta = Q*XTy;
    Sigma_beta = Q/E_ie_2_kappa;
    
    // Update q(kappa)
    E_ie_2_kappa = E_lnig(n + p, mu_kappa, sigma_2_kappa,
                          sum(pow(y - X*mu_beta, 2.0)) + pow(lambda, 2.0)*sum(E_a%pow(mu_beta, 2.0)) + trace(Q_inv*Sigma_beta),
                          "1/exp(2*x)", n_grid);
    
    // Update q(a)
    E_a = sqrt(1/(E_ie_2_kappa*(pow(mu_beta, 2.0) + diagvec(Sigma_beta))))/lambda;
    
    // Check for convergence
    if (i == 0) {
      // Base deltas
      bd_mu = norm(mu_beta - mu_beta_old, 2);
      bd_Sigma = norm(Sigma_beta - Sigma_beta_old, "fro");
      bd_a = norm(E_a - E_a_old, 2);
      
      if (verbose) {
        Rcout << "Delta for mu: " << bd_mu << '\n';
        Rcout << "Delta for Sigma: " << bd_Sigma << '\n';
        Rcout << "Delta for a: " << bd_a << '\n';
      }
      
      continue;
    }
    
    double d_mu = norm(mu_beta - mu_beta_old, 2);
    double d_Sigma = norm(Sigma_beta - Sigma_beta_old, "fro");
    double d_a = norm(E_a - E_a_old, 2);
    
    if (verbose) {
      Rcout << "Delta for mu: " << d_mu << '\n';
      Rcout << "Delta for Sigma: " << d_Sigma << '\n';
      Rcout << "Delta for a: " << d_a << '\n';
    }
    
    if (d_mu < thresh*bd_mu && d_Sigma < thresh*bd_Sigma && d_a < thresh*bd_a && i >= min_iter - 1) {
      converged = true;
      if (verbose) {
        Rcout << "MFVB has converged; stopping MFVB\n";
      }
      break;
    }
  }
  
  if (!converged) {
    stop("MFVB failed to converge");
  }
  
  // Return parameters
  double mu_kappa_q = E_lnig(n + p, mu_kappa, sigma_2_kappa,
                             sum(pow(y - X*mu_beta, 2.0)) + pow(lambda, 2.0)*sum(E_a%pow(mu_beta, 2.0)) + trace(Q_inv*Sigma_beta),
                             "x", n_grid);
  double sigma_2_kappa_q = E_lnig(n + p, mu_kappa, sigma_2_kappa,
                                  sum(pow(y - X*mu_beta, 2.0)) + pow(lambda, 2.0)*sum(E_a%pow(mu_beta, 2.0)) + trace(Q_inv*Sigma_beta),
                                  "x^2", n_grid) - pow(mu_kappa_q, 2.0);
  
  vec mu_theta = zeros(p + 1);
  mu_theta.subvec(0, p - 1) = mu_beta;
  mu_theta(p) = mu_kappa_q;
  
  mat Sigma_theta = zeros(p + 1, p + 1);
  Sigma_theta.submat(0, 0, p - 1, p - 1) = Sigma_beta;
  Sigma_theta(p, p) = sigma_2_kappa_q;
  
  return List::create(_["mu"] = mu_theta, _["Sigma"] = Sigma_theta);
}
