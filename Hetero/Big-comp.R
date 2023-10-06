# Big data EP comparison for heteroscedastic linear regression

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
method <- args[1]
seed <- as.numeric(args[2])
n_grid <- as.numeric(args[3])
min_passes <- as.numeric(args[4])
thresh <- as.numeric(args[5])
set.seed(seed)

load("Hetero/Data/Big/Big.RData")
n <- nrow(X.1)
p.1 <- ncol(X.1)
p.2 <- ncol(X.2)

mu.theta <- rep(0, p.1 + p.2)
Sigma.theta <- diag(c(rep(sigma.2.beta.1, p.1), rep(sigma.2.beta.2, p.2)))

if (method == "ep") {
  start.time <- proc.time()
  
  ep.res <- ep(X.1, X.2, y, Sigma.theta, mu.theta,
               eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
               min_passes = min_passes, max_passes = 200, thresh = thresh, n_grid = n_grid, verbose = T)
  
  total.time <- proc.time() - start.time
  print(total.time["elapsed"])
} else if (method == "ep-2d") {
  start.time <- proc.time()
  
  ep.2d.res <- ep_2d(X.1, X.2, y, Sigma.theta, mu.theta,
                     eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
                     min_passes = min_passes, max_passes = 200, thresh = thresh, n_grid = n_grid, verbose = T)
  
  total.time <- proc.time() - start.time
  print(total.time["elapsed"])
} else {
  stop("method must be in one of ep, or ep-2d")
}
