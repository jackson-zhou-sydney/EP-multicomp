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
  double E;
public:
  Neg_Expnt(double A_, double B_, double C_, double D_, double E_) : A(A_), B(B_), C(C_), D(D_), E(E_) {}
  
  double f_grad(Constvec& x, Refvec grad) {
    double x_dbl = x(0);
    double f = A/exp(2.0*x_dbl) - B/exp(x_dbl) + C*x_dbl + pow(x_dbl - D, 2.0)/(2.0*E);
    grad(0) = -2.0*A/exp(2.0*x_dbl) + B/exp(x_dbl) + C + (x_dbl - D)/E;
    return f;
  }
};

double trap_unif(double delta_x, rowvec y) {
  // Trapezoidal integration with uniform intervals
  y(0) = 0.5*y(0);
  y.back() = 0.5*y.back();
  return delta_x*sum(y);
}

double E_lnig(double A, double B, double C, double D, double E, String fun, int n_grid) {
  // Expectation of log-normal and inverse gamma product
  Neg_Expnt neg_expnt(A, B, C, D, E);
  
  VectorXd x_est = VectorXd::Zero(1);
  double fopt;
  int status = optim_lbfgs(neg_expnt, x_est, fopt, 1000);
  if (status < 0) {
    stop("Failed to converge");
  }
  
  double mu = x_est(0);
  double sigma = 1.0/sqrt((4.0*A*exp(-2.0*mu) - B*exp(-mu) + 1.0/E));
  
  double lb = mu - 5.0*sigma;
  double ub = mu + 5.0*sigma;
  
  vec x_values = linspace(lb, ub, n_grid);
  double delta_x = x_values(1) - x_values(0);
  mat y_matrix = zeros(2, n_grid);
  
  if (fun == "x") {
    for (int i = 0; i < n_grid; ++i) {
      double x = x_values(i);
      y_matrix(0, i) = exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) + fopt);
      y_matrix(1, i) = x*y_matrix(0, i);
    }
  } else if (fun == "x^2") {
    for (int i = 0; i < n_grid; ++i) {
      double x = x_values(i);
      y_matrix(0, i) = exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) + fopt);
      y_matrix(1, i) = pow(x, 2.0)*y_matrix(0, i);
    }
  } else if (fun == "1/exp(x)") {
    for (int i = 0; i < n_grid; ++i) {
      double x = x_values(i);
      y_matrix(0, i) = exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) + fopt);
      y_matrix(1, i) = exp(-x)*y_matrix(0, i);
    }
  } else if (fun == "1/exp(2*x)") {
    for (int i = 0; i < n_grid; ++i) {
      double x = x_values(i);
      y_matrix(0, i) = exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) + fopt);
      y_matrix(1, i) = exp(-2*x)*y_matrix(0, i);
    }
  } else {
    stop("fun must be one of: x, x^2, 1/exp(x), or 1/exp(2*x)");
  }
  
  return trap_unif(delta_x, y_matrix.row(1))/trap_unif(delta_x, y_matrix.row(0));
}

// [[Rcpp::export]]
List mfvb(mat X, vec y, mat Sigma_beta_p, vec mu_beta_p, double sigma_2_kappa, double mu_kappa,
          double tau, int min_iter, int max_iter, double thresh, int n_grid, bool verbose) {
  // MFVB for Bayesian quantile regression
  int n = X.n_rows;
  int p = X.n_cols;
  
  double bd_mu;
  double bd_Sigma;
  double bd_a;
  
  mat Sigma_beta_p_inv = inv_sympd(Sigma_beta_p);
  vec X_colsums = sum(X, 0).t();
  
  // Initialisations
  vec mu_beta = zeros(p);
  mat Sigma_beta = eye(p, p);
  double E_ie_kappa = 1;
  double E_ie_2_kappa = 1;
  vec E_a = ones(n);
  
  vec E_y_X_beta_2 = pow(y - X*mu_beta, 2.0) + diagvec(X*Sigma_beta*X.t());
  mat A = diagmat(E_a);
  
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
    mat Q_inv = E_ie_2_kappa*tau*(1.0 - tau)*X.t()*A*X + Sigma_beta_p_inv;
    mat Q = inv(Q_inv);
    mu_beta = Q*(E_ie_2_kappa*tau*(1.0 - tau)*X.t()*A*y - E_ie_kappa*(0.5 - tau)*X_colsums + Sigma_beta_p_inv*mu_beta_p);
    Sigma_beta = Q;
    E_y_X_beta_2 = pow(y - X*mu_beta, 2.0) + diagvec(X*Sigma_beta*X.t());
    
    // Update q(kappa)
    E_ie_kappa = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                        (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                        "1/exp(x)", n_grid);
    E_ie_2_kappa = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                          (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                          "1/exp(2*x)", n_grid);
    
    // Update q(a)
    E_a = 0.5/(tau*(1 - tau))*pow(E_ie_2_kappa*E_y_X_beta_2, -0.5);
    A = diagmat(E_a);
    
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
      if (verbose) {
        Rcout << "MFVB has converged; stopping MFVB\n";
      }
      break;
    }
  }
  
  // Return parameters
  double mu_kappa_q = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                             (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                             "x", n_grid);
  double sigma_2_kappa_q = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                                  (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                                  "x^2", n_grid) - pow(mu_kappa_q, 2.0);
  
  vec mu_theta = zeros(p + 1);
  mu_theta.subvec(0, p - 1) = mu_beta;
  mu_theta(p) = mu_kappa_q;
  
  mat Sigma_theta = zeros(p + 1, p + 1);
  Sigma_theta.submat(0, 0, p - 1, p - 1) = Sigma_beta;
  Sigma_theta(p, p) = sigma_2_kappa_q;
  
  return List::create(_["mu"] = mu_theta, _["Sigma"] = Sigma_theta);
}
