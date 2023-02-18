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

class MFVB_q: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
public:
  MFVB_q(double A_, double B_, double C_, double D_, double E_, double F_) : A(A_), B(B_), C(C_), D(D_), E(E_), F(F_) {}
  
  double operator()(const double& x) const {
    return exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) - F);
  }
};

class MFVB_r_x: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
public:
  MFVB_r_x(double A_, double B_, double C_, double D_, double E_, double F_) : A(A_), B(B_), C(C_), D(D_), E(E_), F(F_) {}
  
  double operator()(const double& x) const {
    return exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) - F)*x;
  }
};

class MFVB_r_x2: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
public:
  MFVB_r_x2(double A_, double B_, double C_, double D_, double E_, double F_) : A(A_), B(B_), C(C_), D(D_), E(E_), F(F_) {}
  
  double operator()(const double& x) const {
    return exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) - F)*pow(x, 2.0);
  }
};

class MFVB_r_exp_x: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
public:
  MFVB_r_exp_x(double A_, double B_, double C_, double D_, double E_, double F_) : A(A_), B(B_), C(C_), D(D_), E(E_), F(F_) {}
  
  double operator()(const double& x) const {
    return exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) - F - x);
  }
};

class MFVB_r_exp_2x: public Func {
private:
  double A;
  double B;
  double C;
  double D;
  double E;
  double F;
public:
  MFVB_r_exp_2x(double A_, double B_, double C_, double D_, double E_, double F_) : A(A_), B(B_), C(C_), D(D_), E(E_), F(F_) {}
  
  double operator()(const double& x) const {
    return exp(-A/exp(2.0*x) + B/exp(x) - C*x - pow(x - D, 2.0)/(2.0*E) - F - 2.0*x);
  }
};

double E_lnig(double A, double B, double C, double D, double E, String fun) {
  // Expectation of log-normal and inverse gamma product
  Neg_Expnt neg_expnt(A, B, C, D, E);
  
  VectorXd x_est = VectorXd::Zero(1);
  double fopt;
  int status = optim_lbfgs(neg_expnt, x_est, fopt, 1000);
  if (status < 0) {
    stop("Failed to converge");
  }
  
  MFVB_q q(A, B, C, D, E, -fopt);
  double err_est_q;
  double err_est_r;
  int err_code_q;
  int err_code_r;
  
  if (fun == "x") {
    MFVB_r_x r(A, B, C, D, E, -fopt);
    return integrate(r, -5.0, INFINITY, err_est_r, err_code_r, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81)/
           integrate(q, -5.0, INFINITY, err_est_q, err_code_q, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81);
  } else if (fun == "x^2") {
    MFVB_r_x2 r(A, B, C, D, E, -fopt);
    return integrate(r, -5.0, INFINITY, err_est_r, err_code_r, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81)/
           integrate(q, -5.0, INFINITY, err_est_q, err_code_q, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81);
  } else if (fun == "1/exp(x)") {
    MFVB_r_exp_x r(A, B, C, D, E, -fopt);
    return integrate(r, -5.0, INFINITY, err_est_r, err_code_r, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81)/
           integrate(q, -5.0, INFINITY, err_est_q, err_code_q, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81);
  } else if (fun == "1/exp(2*x)") {
    MFVB_r_exp_2x r(A, B, C, D, E, -fopt);
    return integrate(r, -5.0, INFINITY, err_est_r, err_code_r, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81)/
           integrate(q, -5.0, INFINITY, err_est_q, err_code_q, 100, 1e-8, 1e-6, Integrator<double>::GaussKronrod81);
  } else {
    stop("fun must be one of: x, x^2, 1/exp(x), or 1/exp(2*x)");
  }
}

// [[Rcpp::export]]
List mfvb_c(mat X, vec y, mat Sigma_beta_p, vec mu_beta_p, double sigma_2_kappa, double mu_kappa,
            double tau, int maxit, double tol) {
  // MFVB for Bayesian quantile regression
  int n = X.n_rows;
  int p = X.n_cols;
  
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
  for (int i = 0; i < maxit; ++i) {
    // Store old values
    vec mu_beta_old = mu_beta;
    
    // Update q(beta)
    mat Q_inv = E_ie_2_kappa*tau*(1.0 - tau)*X.t()*A*X + Sigma_beta_p_inv;
    mat Q = inv(Q_inv);
    mu_beta = Q*(E_ie_2_kappa*tau*(1.0 - tau)*X.t()*A*y - E_ie_kappa*(0.5 - tau)*X_colsums + Sigma_beta_p_inv*mu_beta_p);
    Sigma_beta = Q;
    E_y_X_beta_2 = pow(y - X*mu_beta, 2.0) + diagvec(X*Sigma_beta*X.t());
    
    // Update q(kappa)
    E_ie_kappa = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                        (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                        "1/exp(x)");
    E_ie_2_kappa = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                          (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                          "1/exp(2*x)");
    
    // Update q(a)
    E_a = 0.5/(tau*(1 - tau))*pow(E_ie_2_kappa*E_y_X_beta_2, -0.5);
    A = diagmat(E_a);
    
    // Check for convergence
    double delta = norm(mu_beta - mu_beta_old);
    if (delta < tol) {
      break;
    }
  }
  
  // Return parameters
  double mu_kappa_q = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                             (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                             "x");
  double sigma_2_kappa_q = E_lnig(0.5*tau*(1.0 - tau)*dot(E_a, E_y_X_beta_2),
                                  (0.5 - tau)*sum(y - X*mu_beta), n, mu_kappa, sigma_2_kappa,
                                  "x^2") - pow(mu_kappa_q, 2.0);
  
  vec mu_theta = zeros(p + 1);
  mu_theta.subvec(0, p - 1) = mu_beta;
  mu_theta(p) = mu_kappa_q;
  
  mat Sigma_theta = zeros(p + 1, p + 1);
  Sigma_theta.submat(0, 0, p - 1, p - 1) = Sigma_beta;
  Sigma_theta(p, p) = sigma_2_kappa_q;
  
  return List::create(_["mu"] = mu_theta, _["Sigma"] = Sigma_theta);
}
