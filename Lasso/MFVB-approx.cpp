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

class MFVB_q: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
public:
  MFVB_q(double A_, double B_, double C_, double D_, double E_) : A(A_), B(B_), C(C_), D(D_), E(E_) {}
  
  double operator()(const double& x) const {
    return exp(-A*x - pow(x - B, 2.0)/(2.0*C) - D/(2.0*exp(2.0*x)) - E);
  }
};

class MFVB_r_x: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
public:
  MFVB_r_x(double A_, double B_, double C_, double D_, double E_) : A(A_), B(B_), C(C_), D(D_), E(E_) {}
  
  double operator()(const double& x) const {
    return exp(-A*x - pow(x - B, 2.0)/(2.0*C) - D/(2.0*exp(2.0*x)) - E)*x;
  }
};

class MFVB_r_x2: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
public:
  MFVB_r_x2(double A_, double B_, double C_, double D_, double E_) : A(A_), B(B_), C(C_), D(D_), E(E_) {}
  
  double operator()(const double& x) const {
    return exp(-A*x - pow(x - B, 2.0)/(2.0*C) - D/(2.0*exp(2.0*x)) - E)*pow(x, 2.0);
  }
};

class MFVB_r_exp_2x: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
public:
  MFVB_r_exp_2x(double A_, double B_, double C_, double D_, double E_) : A(A_), B(B_), C(C_), D(D_), E(E_) {}
  
  double operator()(const double& x) const {
    return exp(-A*x - pow(x - B, 2.0)/(2.0*C) - D/(2.0*exp(2.0*x)) - E - 2*x);
  }
};

double E_lnig_x(double A, double B, double C, double D, String fun) {
  // Expectation of log-normal and inverse gamma product
  Neg_Expnt neg_expnt(A, B, C, D);
  
  VectorXd x_est = VectorXd::Zero(1);
  double fopt;
  int status = optim_lbfgs(neg_expnt, x_est, fopt, 1000);
  if (status < 0) {
    stop("Failed to converge");
  }
  
  MFVB_q q(A, B, C, D, -fopt);
  double err_est_q;
  double err_est_r;
  int err_code_q;
  int err_code_r;
  
  if (fun == "x") {
    MFVB_r_x r(A, B, C, D, -fopt);
    return integrate(r, -INFINITY, INFINITY, err_est_r, err_code_r)/
           integrate(q, -INFINITY, INFINITY, err_est_q, err_code_q);
  } else if (fun == "x^2") {
    MFVB_r_x2 r(A, B, C, D, -fopt);
    return integrate(r, -INFINITY, INFINITY, err_est_r, err_code_r)/
           integrate(q, -INFINITY, INFINITY, err_est_q, err_code_q);
  } else if (fun == "1/exp(2*x)") {
    MFVB_r_exp_2x r(A, B, C, D, -fopt);
    return integrate(r, -INFINITY, INFINITY, err_est_r, err_code_r)/
           integrate(q, -INFINITY, INFINITY, err_est_q, err_code_q);
  } else {
    stop("fun must be one of: x, x^2, or 1/exp(2*x)");
  }
}

// [[Rcpp::export]]
List mfvb_c(mat X, vec y, double sigma_2_kappa, double mu_kappa,
            double lambda, int maxit, double tol) {
  // MFVB for Bayesian lasso linear regression
  int n = X.n_rows;
  int p = X.n_cols;
  
  mat XTX = X.t()*X;
  vec XTy = X.t()*y;
  
  // Initialisations
  vec mu_beta = zeros(p);
  mat Sigma_beta = eye(p, p);
  double E_ie_2_kappa = 1;
  vec E_a = ones(p);
  mat Q;
  mat Q_inv;
  
  // Main MFVB loop
  for (int i = 0; i < maxit; ++i) {
    // Store old values
    vec mu_beta_old = mu_beta;
    
    // Update q(beta)
    Q_inv = XTX + pow(lambda, 2.0)*diagmat(E_a);
    Q = inv(Q_inv);
    mu_beta = Q*XTy;
    Sigma_beta = Q/E_ie_2_kappa;
    
    // Update q(kappa)
    E_ie_2_kappa = E_lnig_x(n + p, mu_kappa, sigma_2_kappa,
                            sum(pow(y - X*mu_beta, 2.0)) + pow(lambda, 2.0)*sum(E_a%pow(mu_beta, 2.0)) + trace(Q_inv*Sigma_beta),
                            "1/exp(2*x)");
    
    // Update q(a)
    E_a = sqrt(1/(E_ie_2_kappa*(pow(mu_beta, 2.0) + diagvec(Sigma_beta))))/lambda;
    
    // Check for convergence
    double delta = norm(mu_beta - mu_beta_old);
    if (delta < tol) {
      break;
    }
  }
  
  // Return parameters
  double mu_kappa_q = E_lnig_x(n + p, mu_kappa, sigma_2_kappa,
                               sum(pow(y - X*mu_beta, 2.0)) + pow(lambda, 2.0)*sum(E_a%pow(mu_beta, 2.0)) + trace(Q_inv*Sigma_beta),
                               "x");
  double sigma_2_kappa_q = E_lnig_x(n + p, mu_kappa, sigma_2_kappa,
                                    sum(pow(y - X*mu_beta, 2.0)) + pow(lambda, 2.0)*sum(E_a%pow(mu_beta, 2.0)) + trace(Q_inv*Sigma_beta),
                                    "x^2") - pow(mu_kappa_q, 2.0);
  
  vec mu_theta = zeros(p + 1);
  mu_theta.subvec(0, p - 1) = mu_beta;
  mu_theta(p) = mu_kappa_q;
  
  mat Sigma_theta = zeros(p + 1, p + 1);
  Sigma_theta.submat(0, 0, p - 1, p - 1) = Sigma_beta;
  Sigma_theta(p, p) = sigma_2_kappa_q;
  
  return List::create(_["mu"] = mu_theta, _["Sigma"] = Sigma_theta);
}
