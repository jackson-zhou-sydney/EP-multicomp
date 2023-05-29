#include <RcppArmadillo.h>
#include <tuple>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double trap_unif_2d(double delta_x, double delta_y, mat Z) {
  // 2D trapezoidal integration with uniform intervals
  Z.row(0) = 0.5*Z.row(0);
  Z.row(Z.n_rows - 1) = 0.5*Z.row(Z.n_rows - 1);
  Z.col(0) = 0.5*Z.col(0);
  Z.col(Z.n_cols - 1) = 0.5*Z.col(Z.n_cols - 1);
  
  return delta_x*delta_y*accu(Z);
}

tuple <vec, mat> h_mom_1_2d(double y, vec mu, mat Sigma, double eta, int n_grid) {
  // Hybrid moments for likelihood sites
  double y_data = y;
  mat Q = inv_sympd(Sigma);
  vec r = Q*mu;
  vec sqrt_dg_Sigma = sqrt(diagvec(Sigma));
  
  vec lb = mu - 5.0*sqrt_dg_Sigma;
  vec ub = mu + 5.0*sqrt_dg_Sigma;
  lb(1) = max(lb(1), -5.0);
  double adjust = 0.5*dot(mu, Q*mu) - dot(mu, r) + eta*(mu(1) + pow(y_data - mu(0), 2.0)/(2.0*exp(2.0*mu(1))));
  
  vec x_values = linspace(lb(0), ub(0), n_grid);
  double delta_x = x_values(1) - x_values(0);
  vec y_values = linspace(lb(1), ub(1), n_grid);
  double delta_y = y_values(1) - y_values(0);
  cube Z_cube = zeros(6, n_grid, n_grid);
  
  for (int i = 0; i < n_grid; ++i) {
    for (int j = 0; j < n_grid; ++j) {
      double x = x_values(i);
      double y = y_values(j);
      vec xy = { x, y };
      
      Z_cube(0, i, j) = exp(-0.5*dot(xy, Q*xy) + dot(xy, r) - eta*(y + pow(y_data - x, 2.0)/(2.0*exp(2.0*y))) + adjust);
      Z_cube(1, i, j) = x*Z_cube(0, i, j);
      Z_cube(2, i, j) = y*Z_cube(0, i, j);
      Z_cube(3, i, j) = pow(x, 2.0)*Z_cube(0, i, j);
      Z_cube(4, i, j) = x*y*Z_cube(0, i, j);
      Z_cube(5, i, j) = pow(y, 2.0)*Z_cube(0, i, j);
    }
  }
  
  double int_0 = trap_unif_2d(delta_x, delta_y, Z_cube.row(0));
  double int_11 = trap_unif_2d(delta_x, delta_y, Z_cube.row(1));
  double int_12 = trap_unif_2d(delta_x, delta_y, Z_cube.row(2));
  double int_211 = trap_unif_2d(delta_x, delta_y, Z_cube.row(3));
  double int_212 = trap_unif_2d(delta_x, delta_y, Z_cube.row(4));
  double int_222 = trap_unif_2d(delta_x, delta_y, Z_cube.row(5));
  
  vec ih1 = { int_11, int_12 };
  
  mat ih2 = { {int_211, int_212}, 
              {int_212, int_222} };
  
  vec mu_h = ih1/int_0;
  mat Sigma_h = ih2/int_0 - mu_h*mu_h.t();
  
  return make_tuple(mu_h, Sigma_h);
}

tuple <vec, mat> h_mom_2_2d(double lambda, vec mu, mat Sigma, double eta, int n_grid) {
  // Hybrid moments for Laplace-based prior sites
  mat Q = inv_sympd(Sigma);
  vec r = Q*mu;
  vec sqrt_dg_Sigma = sqrt(diagvec(Sigma));
  
  vec lb = mu - 5.0*sqrt_dg_Sigma;
  vec ub = mu + 5.0*sqrt_dg_Sigma;
  lb(1) = max(lb(1), -10.0);
  double adjust = 0.5*dot(mu, Q*mu) - dot(mu, r) + eta*(mu(1) + lambda*abs(mu(0))/exp(mu(1)));
  
  vec x_values = linspace(lb(0), ub(0), n_grid);
  double delta_x = x_values(1) - x_values(0);
  vec y_values = linspace(lb(1), ub(1), n_grid);
  double delta_y = y_values(1) - y_values(0);
  cube Z_cube = zeros(6, n_grid, n_grid);
  
  for (int i = 0; i < n_grid; ++i) {
    for (int j = 0; j < n_grid; ++j) {
      double x = x_values(i);
      double y = y_values(j);
      vec xy = { x, y };
      
      Z_cube(0, i, j) = exp(-0.5*dot(xy, Q*xy) + dot(xy, r) - eta*(y + lambda*abs(x)/exp(y)) + adjust);
      Z_cube(1, i, j) = x*Z_cube(0, i, j);
      Z_cube(2, i, j) = y*Z_cube(0, i, j);
      Z_cube(3, i, j) = pow(x, 2.0)*Z_cube(0, i, j);
      Z_cube(4, i, j) = x*y*Z_cube(0, i, j);
      Z_cube(5, i, j) = pow(y, 2.0)*Z_cube(0, i, j);
    }
  }
  
  double int_0 = trap_unif_2d(delta_x, delta_y, Z_cube.row(0));
  double int_11 = trap_unif_2d(delta_x, delta_y, Z_cube.row(1));
  double int_12 = trap_unif_2d(delta_x, delta_y, Z_cube.row(2));
  double int_211 = trap_unif_2d(delta_x, delta_y, Z_cube.row(3));
  double int_212 = trap_unif_2d(delta_x, delta_y, Z_cube.row(4));
  double int_222 = trap_unif_2d(delta_x, delta_y, Z_cube.row(5));
  
  vec ih1 = { int_11, int_12 };
  
  mat ih2 = { {int_211, int_212}, 
              {int_212, int_222} };
  
  vec mu_h = ih1/int_0;
  mat Sigma_h = ih2/int_0 - mu_h*mu_h.t();
  
  return make_tuple(mu_h, Sigma_h);
}

// [[Rcpp::export]]
List ep_2d(mat X, vec y, double sigma_2_kappa, double mu_kappa,
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
          h_mom_res = h_mom_1_2d(y(k), mu_m_star, Sigma_m_star, eta, n_grid);
        } else {
          h_mom_res = h_mom_2_2d(lambda, mu_m_star, Sigma_m_star, eta, n_grid);
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
