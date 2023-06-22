#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <tuple>

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

List laplace(arma::mat X_1, arma::mat X_2, arma::vec y, 
             arma::mat Sigma_theta, arma::vec mu_theta, 
             arma::vec init, int max_iter) {
  // Laplace approximation for Bayesian heteroscedastic linear regression
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

arma::vec vech(arma::mat A) {
  // Convert lower triangular part of matrix to vector
  // A is required to be square
  int d = A.n_rows;
  arma::vec x = arma::zeros(0.5*d*(d + 1));
  for (int i = 0; i < d; ++i) {
    x.subvec(i*(d - 0.5*(i - 1)), (i + 1)*(d - 0.5*i) - 1) = A.submat(i, i, d - 1, i);
  }
  return x;
}

arma::mat vechinv(arma::vec x) {
  // Convert vector to lower triangular matrix
  // The length of x is required to be triangular
  int d = 0.5*(sqrt(8*x.n_elem + 1) - 1);
  arma::mat A = arma::zeros(d, d);
  for (int i = 0; i < d; ++i) {
    A.submat(i, i, d - 1, i) = x.subvec(i*(d - 0.5*(i - 1)), (i + 1)*(d - 0.5*i) - 1);
  }
  return A;
}

double h_lambda(arma::vec theta, arma::mat X_1, arma::mat X_2, arma::vec y, arma::mat Sigma_theta, arma::vec mu_theta, arma::vec lambda) {
  // h_lambda function (up to a constant)
  int n = X_1.n_rows;
  int p_1 = X_1.n_cols;
  int p_2 = X_2.n_cols;
  int d = X_1.n_cols + X_2.n_cols;
  
  arma::vec beta_1 = theta.subvec(0, p_1 - 1);
  arma::vec beta_2 = theta.subvec(p_1, theta.n_elem - 1);
  
  arma::vec A = y - X_1*beta_1;
  arma::vec B = A/exp(2*X_2*beta_2);
  
  arma::vec mu = lambda.subvec(0, d - 1);
  arma::mat L = vechinv(lambda.subvec(d, lambda.n_elem - 1));
  
  return -accu(X_2*beta_2) - 0.5*dot(A, B) - 0.5*dot(theta - mu_theta, inv_sympd(Sigma_theta)*(theta - mu_theta)) + 0.5*dot(theta - mu, inv(L.t())*inv(L)*(theta - mu));
}

arma::vec h_lambda_grad(arma::vec theta, arma::mat X_1, arma::mat X_2, arma::vec y, arma::mat Sigma_theta, arma::vec mu_theta, arma::vec lambda) {
  // Gradient of h_lambda with respect to theta
  int n = X_1.n_rows;
  int p_1 = X_1.n_cols;
  int p_2 = X_2.n_cols;
  int d = X_1.n_cols + X_2.n_cols;
  
  arma::vec beta_1 = theta.subvec(0, p_1 - 1);
  arma::vec beta_2 = theta.subvec(p_1, theta.n_elem - 1);
  
  arma::vec h_grad = arma::zeros(d);
  arma::vec A = y - X_1*beta_1;
  arma::vec B = A/exp(2*X_2*beta_2);
  arma::mat C = X_2;
  for (int i = 0; i < p_2; ++i) {
    C.col(i) = C.col(i)%B;
  }
  h_grad.subvec(0, p_1 - 1) = X_1.t()*B;
  h_grad.subvec(p_1, h_grad.n_elem - 1) = C.t()*A - sum(X_2.t(), 1);
  h_grad = h_grad - inv_sympd(Sigma_theta)*(theta - mu_theta);
  
  arma::vec mu = lambda.subvec(0, d - 1);
  arma::mat L = vechinv(lambda.subvec(d, lambda.n_elem - 1));
  
  return h_grad + inv(L.t())*inv(L)*(theta - mu);
}

// [[Rcpp::export]]
List gvb(arma::mat X_1, arma::mat X_2, arma::vec y, arma::mat Sigma_theta, arma::vec mu_theta,
         int S, arma::vec weights, double eps_0, int tau, arma::vec init, int laplace_max_iter,
         int min_iter, int max_iter, double thresh, bool verbose) {
  // GVB for Bayesian heteroscedastic linear regression
  int n = X_1.n_rows;
  int p_1 = X_1.n_cols;
  int p_2 = X_2.n_cols;
  int d = X_1.n_cols + X_2.n_cols;
  
  double bdelta;
  
  // Initialisation
  List laplace_res = laplace(X_1, X_2, y, Sigma_theta, mu_theta, init, laplace_max_iter);
  arma::mat Sigma_init = laplace_res["Sigma"];
  arma::vec mu_init = laplace_res["mu"];
  
  if (!Sigma_init.is_sympd()) {
    arma::cx_vec eigval_comp;
    arma::cx_mat eigvec_comp;
    eig_gen(eigval_comp, eigvec_comp, Sigma_init);
    arma::vec eigval = real(eigval_comp);
    arma::mat eigvec = real(eigvec_comp);
    
    double min_eigval = abs(eigval).min();
    for (int i = 0; i < eigval.n_elem; ++i) {
      if (eigval(i) <= 0) {
        eigval(i) = min_eigval;
      }
    }
    
    Sigma_init = eigvec*diagmat(eigval)*eigvec.t();
  }
  
  arma::vec lambda = arma::zeros(d + 0.5*d*(d + 1));
  lambda.subvec(d, lambda.n_elem - 1) = vech(chol(Sigma_init, "lower"));
  lambda.subvec(0, d - 1) = mu_init;
  arma::mat L = vechinv(lambda.subvec(d, lambda.n_elem - 1));
  arma::vec mu = lambda.subvec(0, d - 1);
  
  arma::mat mvn_samples = mvnrnd(arma::zeros(d), arma::eye(d, d), S);
  arma::vec g = arma::zeros(d + 0.5*d*(d + 1));
  for (int s = 0; s < S; ++s) {
    arma::vec h_lambda_grad_theta = h_lambda_grad(mu + L*mvn_samples.col(s), X_1, X_2, y, Sigma_theta, mu_theta, lambda);
    g.subvec(0, d - 1) += h_lambda_grad_theta;
    g.subvec(d, g.n_elem - 1) += vech(h_lambda_grad_theta*mvn_samples.col(s).t());
  }
  g /= S;
  arma::vec v = pow(g, 2);
  arma::vec g_bar = g;
  arma::vec v_bar = v;
  
  // Main GVB loop
  for (int t = 0; t < max_iter; ++t) {
    if (verbose) {
      Rcout << "---- Current iteration: " << t << " ----\n";
    }
    
    arma::mat mvn_samples = mvnrnd(arma::zeros(d), arma::eye(d, d), S);
    arma::vec g = arma::zeros(d + 0.5*d*(d + 1));
    for (int s = 0; s < S; ++s) {
      arma::vec h_lambda_grad_theta = h_lambda_grad(mu + L*mvn_samples.col(s), X_1, X_2, y, Sigma_theta, mu_theta, lambda);
      g.subvec(0, d - 1) += h_lambda_grad_theta;
      g.subvec(d, g.n_elem - 1) += vech(h_lambda_grad_theta*mvn_samples.col(s).t());
    }
    g /= S;
    arma::vec v = pow(g, 2);
    g_bar = weights(0)*g_bar + (1 - weights(0))*g;
    v_bar = weights(1)*v_bar + (1 - weights(1))*v;
    
    double alpha;
    if (t < tau) {
      alpha = eps_0;
    } else {
      alpha = eps_0*(tau/(t + 1.0));
    }
    
    lambda += alpha*g_bar/sqrt(v_bar);
    L = vechinv(lambda.subvec(d, lambda.n_elem - 1));
    mu = lambda.subvec(0, d - 1);
    
    double LB = 0.0;
    for (int s = 0; s < S; ++s) {
      LB += h_lambda(mu + L*mvn_samples.col(s), X_1, X_2, y, Sigma_theta, mu_theta, lambda);
    }
    LB /= S;
    
    // Check for convergence
    if (t == 0) {
      // Base delta
      bdelta = norm(alpha*g_bar/sqrt(v_bar), 2);
      
      if (verbose) {
        Rcout << "Delta: " << bdelta << '\n';
        Rcout << "LB: " << LB << '\n';
      }
      
      continue;
    }
    
    double delta = norm(alpha*g_bar/sqrt(v_bar), 2);
    
    if (verbose) {
      Rcout << "Delta: " << delta << '\n';
      Rcout << "LB: " << LB << '\n';
    }
    
    if (delta < thresh*bdelta && t >= min_iter - 1) {
      if (verbose) {
        Rcout << "GVB has converged; stopping GVB\n";
      }
      break;
    }
  }
  
  // Return parameters
  return List::create(_["mu"] = mu, _["Sigma"] = L*L.t());
}
