
# Big data EP comparison for quantile regression

source("General-auxiliaries.R")
source("Quantile/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
method <- args[1]
seed <- as.numeric(args[2])
n_grid <- as.numeric(args[3])
min_passes <- as.numeric(args[4])
thresh <- as.numeric(args[5])
set.seed(seed)

bench.time.df <- data.frame(seed = integer(),
                            bench = integer(),
                            method = character(),
                            n_grid = double(),
                            min_passes = double(),
                            thresh = double(),
                            time = double())

load("Quantile/Data/Big/Big.RData")
n <- nrow(X)
p <- ncol(X)

mu.theta <- rep(0, p + 1)
Sigma.theta <- diag(c(rep(sigma.2.beta, p), sigma.2.kappa))
mu.beta <- mu.theta[1:p]
Sigma.beta <- Sigma.theta[1:p, 1:p]
mu.kappa <- mu.theta[p + 1]
sigma.2.kappa <- Sigma.theta[p + 1, p + 1]

if (method == "ep") {
  start.time <- proc.time()
  
  ep.res <- ep(X, y, Sigma.theta, mu.theta,
               tau, eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
               min_passes = min_passes, max_passes = 200, thresh = thresh, n_grid = n_grid, verbose = T)
  
  total.time <- proc.time() - start.time
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "ep",
                                             n_grid = n_grid,
                                             min_passes = min_passes,
                                             thresh = thresh,
                                             time = total.time["elapsed"])
} else if (method == "ep-2d") {
  start.time <- proc.time()
  
  ep.2d.res <- ep_2d(X, y, Sigma.theta, mu.theta,
                     tau, eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
                     min_passes = min_passes, max_passes = 200, thresh = thresh, n_grid = n_grid, verbose = T)
  
  total.time <- proc.time() - start.time
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "ep-2d",
                                             n_grid = n_grid,
                                             min_passes = min_passes,
                                             thresh = thresh,
                                             time = total.time["elapsed"])
} else {
  stop("method must be in one of ep, or ep-2d")
}

save(bench.time.df,
     file = paste0("Quantile/Results/Big-comp-results-", toupper(method), "-", str_pad(seed, 2, pad = "0"), "-", n_grid, "-", min_passes, "-", thresh,  ".RData"))
