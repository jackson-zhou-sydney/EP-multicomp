#include <RcppArmadillo.h>
#include <tuple>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

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

tuple <vec, mat> h_mom_1(double y, vec mu, mat Sigma, double eta, int n_grid) {
  // Hybrid moments for likelihood sites
  double mu_1 = mu(0);
  double mu_2 = mu(1);
  
  mat Q = inv_sympd(Sigma);
  double Q_11 = Q(0, 0);
  double Q_12 = Q(0, 1);
  double Q_22 = Q(1, 1);
  
  double lb = mu_2 - 5.0*sqrt(Sigma(1, 1));
  double ub = mu_2 + 5.0*sqrt(Sigma(1, 1));
  double adjust = -eta*(2.0*mu_2 + pow(y - mu_1, 2)/exp(2.0*mu_2));
  
  vec x_values = linspace(lb, ub, n_grid);
  double delta_x = x_values(1) - x_values(0);
  mat y_matrix = zeros(6, n_grid);
  
  for (int i = 0; i < n_grid; ++i) {
    double x = x_values(i);
    double a = Q_11 + eta/exp(2.0*x);
    double b = 2.0*(Q_12*(x - mu_2) - Q_11*mu_1) - eta*2.0*y/exp(2.0*x);
    double c = Q_11*pow(mu_1, 2.0) + 2.0*Q_12*mu_1*(mu_2 - x) + Q_22*pow(x - mu_2, 2.0) + eta*(2.0*x + pow(y, 2.0)/exp(2.0*x)) + adjust;
    
    y_matrix(0, i) = gi_0(a, b, c);
    y_matrix(1, i) = gi_1(a, b, c);
    y_matrix(2, i) = x*y_matrix(0, i);
    y_matrix(3, i) = gi_2(a, b, c);
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
  
  vec mu_h = ih1/int_0;
  mat Sigma_h = ih2/int_0 - mu_h*mu_h.t();
  
  return make_tuple(mu_h, Sigma_h);
}

tuple <vec, mat> h_mom_2(double lambda, vec mu, mat Sigma, double eta, int n_grid) {
  // Hybrid moments for Laplace-based prior sites
  double mu_1 = mu(0);
  double mu_2 = mu(1);
  
  mat Q = inv_sympd(Sigma);
  double Q_11 = Q(0, 0);
  double Q_12 = Q(0, 1);
  double Q_22 = Q(1, 1);
  
  double lb = mu_2 - 5.0*sqrt(Sigma(1, 1));
  double ub = mu_2 + 5.0*sqrt(Sigma(1, 1));
  double adjust = -eta*2.0*mu_2;
  
  vec x_values = linspace(lb, ub, n_grid);
  double delta_x = x_values(1) - x_values(0);
  mat y_matrix = zeros(6, n_grid);
  
  for (int i = 0; i < n_grid; ++i) {
    double x = x_values(i);
    double a = Q_11;
    double b_l = 2.0*(Q_12*(x - mu_2) - Q_11*mu_1) - eta*2.0*lambda/exp(x);
    double b_u = 2.0*(Q_12*(x - mu_2) - Q_11*mu_1) + eta*2.0*lambda/exp(x);
    double c = Q_11*pow(mu_1, 2.0) + 2.0*Q_12*mu_1*(mu_2 - x) + Q_22*pow(x - mu_2, 2.0) + eta*2.0*x + adjust;
    
    y_matrix(0, i) = tgi_m_0(a, b_l, c, 0.0) + tgi_p_0(a, b_u, c, 0.0);
    y_matrix(1, i) = tgi_m_1(a, b_l, c, 0.0) + tgi_p_1(a, b_u, c, 0.0);
    y_matrix(2, i) = x*y_matrix(0, i);
    y_matrix(3, i) = tgi_m_2(a, b_l, c, 0.0) + tgi_p_2(a, b_u, c, 0.0);
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
  
  vec mu_h = ih1/int_0;
  mat Sigma_h = ih2/int_0 - mu_h*mu_h.t();
  
  return make_tuple(mu_h, Sigma_h);
}

// [[Rcpp::export]]
List ep(mat X, vec y, double sigma_2_kappa, double mu_kappa,
        double lambda, double eta, double alpha, mat Q_star_init, vec r_star_init,
        double min_passes, double max_passes, double thresh, int n_grid, bool verbose) {
  // Dampened power EP for Bayesian lasso linear regression
  int n = X.n_rows;
  int p = X.n_cols;
  
  double bmd_Q;
  double bmd_r;
  
  // Parameter initialisation
  cube Q_star_values = zeros(2, 2, n + p);
  mat r_star_values = zeros(2, n + p);
  
  mat Q_p = zeros(p + 1, p + 1);
  Q_p(p, p) = 1/sigma_2_kappa;
  vec r_p = zeros(p + 1);
  r_p(p) = mu_kappa/sigma_2_kappa;
  
  mat Q_dot = Q_p;
  vec r_dot = r_p;
  
  mat sym_Q_star_init = symmatu(Q_star_init);
  
  for (int k = 0; k < n + p; ++k) {
    Q_star_values.slice(k) = sym_Q_star_init;
    r_star_values.col(k) = r_star_init;
    
    mat A = zeros(p + 1, 2);
    if (k < n) {
      A.submat(0, 0, p - 1, 0) = X.row(k).t();
      A(p, 1) = 1.0;
    } else {
      A(k - n, 0) = 1.0;
      A(p, 1) = 1.0;
    }
    
    Q_dot += A*Q_star_values.slice(k)*A.t();
    r_dot += A*r_star_values.col(k);
  }
  
  mat Sigma_dot = symmatu(inv_sympd(Q_dot));
  vec mu_dot = Sigma_dot*r_dot;
  
  // Main EP loop
  for (int pass = 0; pass < max_passes; ++pass) {
    // Delta initialisation
    vec deltas_Q = zeros(n + p);
    vec deltas_r = zeros(n + p);
    
    if (verbose) {
      Rcout << "---- Current pass: " << pass << " ----\n";
    }
    
    Q_dot = Q_p;
    r_dot = r_p;
    
    #pragma omp parallel
    {
      mat Q_dot_local = zeros(p + 1, p + 1);
      vec r_dot_local = zeros(p + 1);
      
      #pragma omp for nowait 
      for (int k = 0; k < n + p; ++k) {
        mat A = zeros(p + 1, 2);
        if (k < n) {
          A.submat(0, 0, p - 1, 0) = X.row(k).t();
          A(p, 1) = 1.0;
        } else {
          A(k - n, 0) = 1.0;
          A(p, 1) = 1.0;
        }
        
        mat Sigma_dot_star = A.t()*Sigma_dot*A;
        vec mu_dot_star = A.t()*mu_dot;
        
        mat Q_dot_star = symmatu(inv_sympd(Sigma_dot_star));
        vec r_dot_star = Q_dot_star*mu_dot_star;
        
        mat Q_m_star = Q_dot_star - eta*Q_star_values.slice(k);
        vec r_m_star = r_dot_star - eta*r_star_values.col(k);
        
        mat Sigma_m_star = symmatu(inv_sympd(Q_m_star));
        vec mu_m_star = Sigma_m_star*r_m_star;
        
        tuple <vec, mat> h_mom_res;
        if (k < n) {
          h_mom_res = h_mom_1(y(k), mu_m_star, Sigma_m_star, eta, n_grid);
        } else {
          h_mom_res = h_mom_2(lambda, mu_m_star, Sigma_m_star, eta, n_grid);
        }
        
        mat Sigma_h_star = get<1>(h_mom_res);
        vec mu_h_star = get<0>(h_mom_res);
        
        mat Q_h_star = symmatu(inv_sympd(Sigma_h_star));
        vec r_h_star = Q_h_star*mu_h_star;
        
        mat Q_star_tilde = (1.0 - alpha)*Q_star_values.slice(k) + (alpha/eta)*(Q_h_star - Q_m_star);
        vec r_star_tilde = (1.0 - alpha)*r_star_values.col(k) + (alpha/eta)*(r_h_star - r_m_star);
        
        mat Q_star_tilde_d = Q_star_tilde - Q_star_values.slice(k);
        vec r_star_tilde_d = r_star_tilde - r_star_values.col(k);
        
        deltas_Q(k) = norm(Q_star_tilde_d, "fro");
        deltas_r(k) = norm(r_star_tilde_d, 2);
        
        Q_star_values.slice(k) = Q_star_tilde;
        r_star_values.col(k) = r_star_tilde;
        
        Q_dot_local += A*Q_star_values.slice(k)*A.t();
        r_dot_local += A*r_star_values.col(k);
      }
      
      #pragma omp critical
      {
        Q_dot += Q_dot_local;
        r_dot += r_dot_local;
      }
    }
    
    Sigma_dot = symmatu(inv_sympd(symmatu(Q_dot)));
    mu_dot = Sigma_dot*r_dot;
    
    if (pass == 0) {
      // Base maximum deltas
      bmd_Q = max(deltas_Q);
      bmd_r = max(deltas_r);
      
      if (verbose) {
        Rcout << "Maximum delta for Q: " << bmd_Q << '\n';
        Rcout << "Maximum delta for r: " << bmd_r << '\n';
      }
      
      continue;
    }
    
    double md_Q = max(deltas_Q);
    double md_r = max(deltas_r);
    
    if (verbose) {
      Rcout << "Maximum delta for Q: " << md_Q << '\n';
      Rcout << "Maximum delta for r: " << md_r << '\n';
    }
    
    if (md_Q < thresh*bmd_Q && md_r < thresh*bmd_r && pass > min_passes) {
      if (verbose) {
        Rcout << "EP has converged; stopping EP\n";
      }
      break;
    }
  }
  
  // Returning in moment parameterisation
  return List::create(_["mu"] = mu_dot, _["Sigma"] = Sigma_dot);
}
