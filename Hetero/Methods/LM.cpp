#include <RcppArmadillo.h>
#include <RcppNumerical.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

using namespace Rcpp;
using namespace Eigen;
using namespace Numer;
using namespace std;

arma::vec cast_arma_vec(VectorXd x) {
  // Convert Eigen vector to Armadillo vector
  return arma::vec(x.data(), x.rows(), false, false);
}

arma::mat cast_arma_mat(MatrixXd A) {
  // Convert Eigen matrix to Armadillo matrix
  return arma::mat(A.data(), A.rows(), A.cols(), false, false);
}

typedef Map<MatrixXd> MapMat;
typedef Map<VectorXd> MapVec;

class HeteroReg: public MFuncGrad {
private:
  const MapMat M_X_1;
  const MapMat M_X_2;
  const MapVec M_y;
  const MapMat M_Sigma_theta;
  const MapVec M_mu_theta;
  int n;
  int p_1;
  int p_2;
public:
  HeteroReg(const MapMat M_X_1_, const MapMat M_X_2_, const MapVec M_y_,
            const MapMat M_Sigma_theta_, const MapVec M_mu_theta_) : 
  M_X_1(M_X_1_), 
  M_X_2(M_X_2_), 
  M_y(M_y_), 
  M_Sigma_theta(M_Sigma_theta_), 
  M_mu_theta(M_mu_theta_) {
    n = M_X_1_.rows();
    p_1 = M_X_1_.cols();
    p_2 = M_X_2_.cols();
  }
  
  double f_grad(Constvec& theta, Refvec grad) {
    // Negative log joint likelihood
    VectorXd beta_1 = VectorXd::Zero(p_1);
    for (int i = 0; i < p_1; ++i) {
      beta_1(i) = theta(i);
    }
    
    VectorXd beta_2 = VectorXd::Zero(p_2);
    for (int i = 0; i < p_2; ++i) {
      beta_2(i) = theta(p_1 + i);
    }
    
    VectorXd yXbXb = (M_y - M_X_1*beta_1).array()/exp((2*M_X_2*beta_2).array());
    VectorXd theta_mu_theta = theta - M_mu_theta;
    MatrixXd Q_theta = M_Sigma_theta.inverse();
    
    double f = 0.5*((2*M_X_2*beta_2).sum() + 
                    (M_y - M_X_1*beta_1).dot(yXbXb) +
                    theta_mu_theta.transpose()*Q_theta*theta_mu_theta);
    
    // Gradient
    VectorXd a_1 = M_X_1.transpose()*yXbXb;
    MatrixXd X_3 = MatrixXd::Zero(n, p_2);
    for (int j = 0; j < p_2; ++j) {
      X_3.col(j) = 2*M_X_2.col(j).array()*yXbXb.array();
    }
    VectorXd a_2 = 0.5*(X_3.transpose()*(M_y - M_X_1*beta_1) - (VectorXd)(2*M_X_2.colwise().sum()));
    
    VectorXd a(p_1 + p_2);
    a << a_1, a_2;
    
    grad = Q_theta*theta - a;
    
    return f;
  }
};

// [[Rcpp::export]]
List lm(arma::mat X_1, arma::mat X_2, arma::vec y, 
        arma::mat Sigma_theta, arma::vec mu_theta, 
        arma::vec init, int max_iter) {
  // Laplace's method for Bayesian heteroscedastic linear regression
  const MapMat M_X_1 = as<MapMat>(wrap(X_1));
  const MapMat M_X_2 = as<MapMat>(wrap(X_2));
  const MapVec M_y = as<MapVec>(wrap(y));
  const MapMat M_Sigma_theta = as<MapMat>(wrap(Sigma_theta));
  const MapVec M_mu_theta = as<MapVec>(wrap(mu_theta));
  const MapVec M_init = as<MapVec>(wrap(init));
  
  HeteroReg nll(M_X_1, M_X_2, M_y, M_Sigma_theta, M_mu_theta);
  
  VectorXd theta_est = M_init;
  double fopt;
  int status = optim_lbfgs(nll, theta_est, fopt, max_iter);
  if (status < 0) {
    stop("Failed to converge");
  }
  
  // Calculate the covariance matrix
  int n = M_X_1.rows();
  int p_1 = M_X_1.cols();
  int p_2 = M_X_2.cols();
  
  VectorXd beta_1 = VectorXd::Zero(p_1);
  for (int i = 0; i < p_1; ++i) {
    beta_1(i) = theta_est(i);
  }
  
  VectorXd beta_2 = VectorXd::Zero(p_2);
  for (int i = 0; i < p_2; ++i) {
    beta_2(i) = theta_est(p_1 + i);
  }
  
  VectorXd yXb = M_y - M_X_1*beta_1;
  VectorXd exp_Xb = exp((2*M_X_2*beta_2).array());
  VectorXd yXbXb = yXb.array()/exp_Xb.array();
  
  MatrixXd X_3 = MatrixXd::Zero(n, p_2);
  for (int j = 0; j < p_2; ++j) {
    X_3.col(j) = 2*M_X_2.col(j).array()*yXbXb.array();
  }
  
  MatrixXd X_4 = MatrixXd::Zero(n, p_1);
  for (int i = 0; i < n; ++i) {
    X_4.row(i) = M_X_1.row(i)/exp_Xb(i);
  }
  
  MatrixXd X_5 = MatrixXd::Zero(n, p_2);
  for (int j = 0; j < p_2; ++j) {
    X_5.col(j) = 2*M_X_2.col(j).array()*yXb.array();
  }
  
  MatrixXd A = M_X_1.transpose()*X_4;
  MatrixXd B = M_X_1.transpose()*X_3;
  MatrixXd D = 0.5*X_5.transpose()*X_3;
  
  MatrixXd AB(p_1, p_1 + p_2);
  MatrixXd CD(p_2, p_1 + p_2);
  MatrixXd ABCD(p_1 + p_2, p_1 + p_2);
  
  AB << A, B;
  CD << B.transpose(), D;
  ABCD << AB, CD;
  
  return List::create(_["mu"] = cast_arma_vec(theta_est), _["Sigma"] = cast_arma_mat((ABCD + M_Sigma_theta.inverse()).inverse()));
}
