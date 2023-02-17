#include <RcppArmadillo.h>
#include <tuple>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double trap_unif(double delta_x, rowvec y) {
  // Trapezoidal integration with uniform intervals
  y(0) = 0.5*y(0);
  y.back() = 0.5*y.back();
  return delta_x*sum(y);
}

double gi_0(double a, double b, double c) {
  // Gaussian integral (0th raw moment)
  return sqrt(M_2PI/a)*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)));
}

double gi_1(double a, double b, double c) {
  // Gaussian integral (1st raw moment)
  return sqrt(M_2PI/a)*(-b/(2.0*a))*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)));
}

double gi_2(double a, double b, double c) {
  // Gaussian integral (2nd raw moment)
  return sqrt(M_2PI/a)*(1/a + pow(b, 2.0)/(4.0*pow(a, 2.0)))*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)));
}

double tgi_m_0(double a, double b, double c, double L) {
  // Truncated Gaussian integral (lower, 0th raw moment)
  return exp(log(sqrt(M_2PI/a)) - 0.5*(c - pow(b, 2.0)/(4.0*a)) + R::pnorm(b/(2.0*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 1, 1));
}

double tgi_m_1(double a, double b, double c, double L) {
  // Truncated Gaussian integral (lower, 1st raw moment)
  return sqrt(M_2PI/a)*(-(b/(2.0*a))*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)) + R::pnorm(b/(2.0*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 1, 1)) - (1/sqrt(a))*exp(-0.5*(c - pow(b, 2.0)/(4*a)) + R::dnorm(b/(2*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 1)));
}

double tgi_m_2(double a, double b, double c, double L) {
  // Truncated Gaussian integral (lower, 2nd raw moment)
  return sqrt(M_2PI/a)*((pow(b, 2.0)/(4.0*a) + 1)*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)) + R::pnorm(b/(2.0*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 1, 1)) + (b/(2.0*sqrt(a)) - sqrt(a)*L)*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)) + R::dnorm(b/(2*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 1)))/a;
}

double tgi_p_0(double a, double b, double c, double L) {
  // Truncated Gaussian integral (upper, 0th raw moment)
  return exp(log(sqrt(M_2PI/a)) - 0.5*(c - pow(b, 2.0)/(4.0*a)) + R::pnorm(b/(2.0*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 0, 1));
}

double tgi_p_1(double a, double b, double c, double L) {
  // Truncated Gaussian integral (upper, 1st raw moment)
  return sqrt(M_2PI/a)*(-(b/(2.0*a))*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)) + R::pnorm(-b/(2.0*sqrt(a)) - sqrt(a)*L, 0.0, 1.0, 1, 1)) + (1/sqrt(a))*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)) + R::dnorm(b/(2.0*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 1)));
}

double tgi_p_2(double a, double b, double c, double L) {
  // Truncated Gaussian integral (upper, 2nd raw moment)
  return sqrt(M_2PI/a)*((pow(b, 2.0)/(4.0*a) + 1)*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)) + R::pnorm(-b/(2.0*sqrt(a)) - sqrt(a)*L, 0.0, 1.0, 1, 1)) - (b/(2.0*sqrt(a)) - sqrt(a)*L)*exp(-0.5*(c - pow(b, 2.0)/(4.0*a)) + R::dnorm(b/(2.0*sqrt(a)) + sqrt(a)*L, 0.0, 1.0, 1)))/a;
}

tuple <double, vec, mat> ihr(double y, double tau, vec mu, mat Sigma, double eta) {
  // Hybrid integrals for likelihood sites
  double mu_1 = mu(0);
  double mu_2 = mu(1);
  
  mat Q = symmatu(inv_sympd(Sigma));
  double Q_11 = Q(0, 0);
  double Q_12 = Q(0, 1);
  double Q_22 = Q(1, 1);
  
  double lb = max(mu_2 - 5.0*sqrt(Sigma(1, 1)), -10.0);
  double ub = mu_2 + 5.0*sqrt(Sigma(1, 1));
  
  vec x_values = linspace(lb, ub, 50);
  double delta_x = x_values(1) - x_values(0);
  mat y_matrix = zeros(6, 50);
  
  for (int i = 0; i < 50; ++i) {
    double x = x_values(i);
    double a = Q_11;
    double b_l = 2.0*(Q_12*(x - mu_2) - Q_11*mu_1) + eta*(-2.0*tau)/exp(x);
    double b_u = 2.0*(Q_12*(x - mu_2) - Q_11*mu_1) + eta*(2.0 - 2.0*tau)/exp(x);
    double c_l = Q_11*pow(mu_1, 2.0) + 2.0*Q_12*mu_1*(mu_2 - x) + Q_22*pow(x - mu_2, 2.0) + eta*(2.0*x + (2.0*tau)*y/exp(x));
    double c_u = Q_11*pow(mu_1, 2.0) + 2.0*Q_12*mu_1*(mu_2 - x) + Q_22*pow(x - mu_2, 2.0) + eta*(2.0*x + (2.0*tau - 2.0)*y/exp(x));
    
    y_matrix(0, i) = tgi_m_0(a, b_l, c_l, y) + tgi_p_0(a, b_u, c_u, y);
    y_matrix(1, i) = tgi_m_1(a, b_l, c_l, y) + tgi_p_1(a, b_u, c_u, y);
    y_matrix(2, i) = x*y_matrix(0, i);
    y_matrix(3, i) = tgi_m_2(a, b_l, c_l, y) + tgi_p_2(a, b_u, c_u, y);
    y_matrix(4, i) = x*y_matrix(1, i);
    y_matrix(5, i) = pow(x, 2.0)*y_matrix(0, i);
  }
  
  double int_0 = trap_unif(delta_x, y_matrix.row(0));
  double int_11 = trap_unif(delta_x, y_matrix.row(1));
  double int_12 = trap_unif(delta_x, y_matrix.row(2));
  double int_211 = trap_unif(delta_x, y_matrix.row(3));
  double int_212 = trap_unif(delta_x, y_matrix.row(4));
  double int_222 = trap_unif(delta_x, y_matrix.row(5));
  
  vec ih1 = { int_11, int_12 };
  
  mat ih2 = { {int_211, int_212}, 
              {int_212, int_222} };
  
  double C = sqrt(det(M_2PI*Sigma));
  
  return make_tuple(int_0/C, ih1/C, ih2/C);
}

// [[Rcpp::export]]
List ep_c(mat X, vec y, mat Sigma_theta, vec mu_theta,
          double tau, double eta, double alpha, mat Q_star_init, vec r_star_init, mat offset,
          double min_passes, double max_passes, double tol, double stop,
          double abs_thresh, double rel_thresh, double delta_limit, double patience) {
  // Dampened power EP for Bayesian quantile regression
  int n = X.n_rows;
  int p = X.n_cols;
  int pat_count = 0;
  
  // Parameter initialisation
  cube Q_star_values = zeros(2, 2, n);
  mat r_star_values = zeros(2, n);
  
  mat Q_p = symmatu(inv_sympd(Sigma_theta));
  vec r_p = Q_p * mu_theta;
  
  mat Q_dot = symmatu(offset) + Q_p;
  vec r_dot = r_p;
  
  mat sym_Q_star_init = symmatu(Q_star_init);
  
  for (int k = 0; k < n; ++k) {
    Q_star_values.slice(k) = sym_Q_star_init;
    r_star_values.col(k) = r_star_init;
    
    mat A = zeros(p + 1, 2);
    A.submat(0, 0, p - 1, 0) = X.row(k).t();
    A(p, 1) = 1.0;
    
    Q_dot += symmatu(A*Q_star_values.slice(k)*A.t());
    r_dot += A*r_star_values.col(k);
  }
  
  mat Sigma_dot = symmatu(inv_sympd(Q_dot));
  vec mu_dot = Sigma_dot*r_dot;
  
  // Delta initialisation
  mat deltas = zeros(max_passes*n, 3);
  double prev_max_delta = INFINITY;
  double prev_med_delta = INFINITY;
  int ind = 0;
  
  // Main EP loop
  for (int pass = 0; pass < max_passes; ++pass) {
    int pass_start_ind = ind;
    
    for (int k : randperm(n)) {
      mat A = zeros(p + 1, 2);
      A.submat(0, 0, p - 1, 0) = X.row(k).t();
      A(p, 1) = 1.0;
      
      mat Sigma_dot_star = symmatu(A.t()*Sigma_dot*A);
      vec mu_dot_star = A.t()*mu_dot;
      
      mat Q_dot_star = symmatu(inv_sympd(Sigma_dot_star));
      vec r_dot_star = Q_dot_star*mu_dot_star;
      
      mat Q_m_star = Q_dot_star - eta*Q_star_values.slice(k);
      vec r_m_star = r_dot_star - eta*r_star_values.col(k);
      
      mat Sigma_m_star = symmatu(inv_sympd(Q_m_star));
      vec mu_m_star = Sigma_m_star*r_m_star;
      
      tuple <double, vec, mat> ihr_res = ihr(y(k), tau, mu_m_star, Sigma_m_star, eta);
      double ih0 = get<0>(ihr_res);
      vec ih1 = get<1>(ihr_res);
      mat ih2 = get<2>(ihr_res);
      
      mat Sigma_h_star = symmatu(ih2/ih0 - ih1*ih1.t()/pow(ih0, 2.0));
      vec mu_h_star = ih1/ih0;
      
      mat Q_h_star = symmatu(inv_sympd(Sigma_h_star));
      vec r_h_star = Q_h_star*mu_h_star;
      
      mat Q_star_tilde = (1.0 - alpha)*Q_star_values.slice(k) + (alpha/eta)*(Q_h_star - Q_m_star);
      vec r_star_tilde = (1.0 - alpha)*r_star_values.col(k) + (alpha/eta)*(r_h_star - r_m_star);
      
      mat Q_star_tilde_d = Q_star_tilde - Q_star_values.slice(k);
      vec r_star_tilde_d = r_star_tilde - r_star_values.col(k);
      
      double delta = max(norm(Q_star_tilde_d, "fro"), norm(r_star_tilde_d, 2));
      
      if (delta > tol*prev_med_delta || delta > delta_limit) {
        continue;
      } else {
        deltas.row(ind) = { (double)pass, (double)k, delta };
        ind += 1;
      }
      
      Q_star_values.slice(k) = Q_star_tilde;
      r_star_values.col(k) = r_star_tilde;
      
      Q_dot += symmatu(A*Q_star_tilde_d*A.t());
      r_dot += A*r_star_tilde_d;
      
      mat Sigma_dot_A = Sigma_dot*A;
      Sigma_dot -= symmatu(Sigma_dot_A*inv(inv(Q_star_tilde_d) + Sigma_dot_star)*Sigma_dot_A.t());
      mu_dot = Sigma_dot*r_dot;
    }
    
    vec pass_deltas = deltas.submat(pass_start_ind, 2, ind - 1, 2);
    double max_delta = max(pass_deltas);
    double med_delta = median(pass_deltas);
    
    if (max_delta < abs_thresh && pass > min_passes) {
      break;
    } else if (max_delta > stop*prev_max_delta) {
      Rcout << "Unstable deltas; stopping EP";
      break;
    }
    
    if (max_delta > rel_thresh*prev_max_delta) {
      pat_count += 1;
    } else {
      pat_count = 0;
    }
    
    if (pat_count == patience) {
      Rcout << "Out of patience; stopping EP";
      break;
    }
    
    prev_max_delta = max_delta;
    prev_med_delta = med_delta;
  }
  
  // Returning in moment parameterisation
  return List::create(_["mu"] = mu_dot, _["Sigma"] = Sigma_dot);
}
