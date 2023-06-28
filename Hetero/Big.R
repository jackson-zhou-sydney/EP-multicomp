
# Big data example for heteroscedastic linear regression

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
method <- args[1]
seed <- as.numeric(args[2])
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Hetero/Methods/MCMC.stan", stanc_options = list("O1"))

temp.env <- new.env()

bench.l1.df <- data.frame(seed = integer(),
                          bench = integer(),
                          method = character(),
                          j = integer(),
                          l1 = double())

bench.mmd.df <- data.frame(seed = integer(),
                           bench = integer(),
                           method = character(),
                           mmd = double())

bench.lppd.df <- data.frame(seed = integer(),
                            bench = integer(),
                            method = character(),
                            lppd = double())

bench.cov.norm.df <- data.frame(seed = integer(),
                                bench = integer(),
                                method = character(),
                                cov_norm = double())

bench.time.df <- data.frame(seed = integer(),
                            bench = integer(),
                            method = character(),
                            time = double())

load("Hetero/Data/Big/Big.RData")
n <- nrow(X.1)
p.1 <- ncol(X.1)
p.2 <- ncol(X.2)

mu.theta <- rep(0, p.1 + p.2)
Sigma.theta <- diag(c(rep(sigma.2.beta.1, p.1), rep(sigma.2.beta.2, p.2)))

load(paste0("Hetero/Data/Big/Big-test-", str_pad(seed, 2, pad = "0"), ".RData"))

load(paste0("Hetero/Results/Big-MCMC-results-", str_pad(seed, 2, pad = "0"), ".Rdata"), envir = temp.env)
mcmc.g.mu <- temp.env$mcmc.g.mu
mcmc.g.Sigma <- temp.env$mcmc.g.Sigma
tail.mcmc.g.samples <- temp.env$tail.mcmc.g.samples
grid.points <- temp.env$grid.points
mcmc.g.values <- temp.env$mcmc.g.values

if (method == "mcmc") {
  load(paste0("Hetero/Results/Big-MCMC-results-", str_pad(seed + 1, 2, pad = "0"), ".Rdata"), envir = temp.env)
  mcmc.mu <- temp.env$mcmc.g.mu
  mcmc.Sigma <- temp.env$mcmc.g.Sigma
  tail.mcmc.samples <- temp.env$tail.mcmc.g.samples
  mcmc.grid.points <- temp.env$grid.points
  mcmc.values <- temp.env$mcmc.g.values
  mcmc.time <- temp.env$mcmc.g.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = 4,
                                           method = "mcmc",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], abs(mcmc.g.values[j, ] - approx(mcmc.grid.points[j, ], mcmc.values[j, ], grid.points[j, ])$y))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = 4,
                                                                 method = "mcmc",
                                                                 mmd = max(kmmd(tail.mcmc.samples, 
                                                                                tail.mcmc.g.samples)@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = 4,
                                                     method = "mcmc",
                                                     cov_norm = norm(mcmc.g.Sigma - mcmc.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "mcmc",
                                             time = mcmc.time)
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "mcmc",
                                             lppd = lppd(X.1.test, X.2.test, y.test, tail.mcmc.samples))
} else if (method == "mcmc-s") {
  load(paste0("Hetero/Results/Big-MCMC-results-", str_pad(seed + 2, 2, pad = "0"), ".Rdata"), envir = temp.env)
  mcmc.s.mu <- temp.env$mcmc.s.mu
  mcmc.s.Sigma <- temp.env$mcmc.s.Sigma
  tail.mcmc.s.samples <- temp.env$tail.mcmc.s.samples
  mcmc.s.grid.points <- temp.env$grid.points
  mcmc.s.values <- temp.env$mcmc.s.values
  mcmc.s.time <- temp.env$mcmc.s.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = 4,
                                           method = "mcmc-s",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], abs(mcmc.g.values[j, ] - approx(mcmc.s.grid.points[j, ], mcmc.s.values[j, ], grid.points[j, ])$y))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = 4,
                                                                 method = "mcmc-s",
                                                                 mmd = max(kmmd(tail.mcmc.s.samples, 
                                                                                tail.mcmc.g.samples)@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = 4,
                                                     method = "mcmc-s",
                                                     cov_norm = norm(mcmc.g.Sigma - mcmc.s.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "mcmc-s",
                                             time = mcmc.s.time)
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "mcmc-s",
                                             lppd = lppd(X.1.test, X.2.test, y.test, tail.mcmc.s.samples))
} else if (method == "ep") {
  start.time <- proc.time()
  
  ep.res <- ep(X.1, X.2, y, Sigma.theta, mu.theta,
               eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
               min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = 4,
                                           method = "ep",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], 
                                                          abs(mcmc.g.values[j, ] - dnorm(grid.points[j, ], ep.mu[j], sqrt(ep.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = 4,
                                                                 method = "ep",
                                                                 mmd = max(kmmd(ep.samples, 
                                                                                tail.mcmc.g.samples)@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = 4,
                                                     method = "ep",
                                                     cov_norm = norm(mcmc.g.Sigma - ep.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "ep",
                                             time = total.time["elapsed"])
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "ep",
                                             lppd = lppd(X.1.test, X.2.test, y.test, ep.samples))
} else if (method == "ep-2d") {
  start.time <- proc.time()
  
  ep.2d.res <- ep_2d(X.1, X.2, y, Sigma.theta, mu.theta,
                     eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
                     min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
  ep.2d.mu <- ep.2d.res$mu
  ep.2d.Sigma <- ep.2d.res$Sigma
  ep.2d.samples <- rmvnorm(eval.size, ep.2d.mu, ep.2d.Sigma)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = 4,
                                           method = "ep-2d",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], 
                                                          abs(mcmc.g.values[j, ] - dnorm(grid.points[j, ], ep.2d.mu[j], sqrt(ep.2d.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = 4,
                                                                 method = "ep-2d",
                                                                 mmd = max(kmmd(ep.2d.samples, 
                                                                                tail.mcmc.g.samples)@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = 4,
                                                     method = "ep-2d",
                                                     cov_norm = norm(mcmc.g.Sigma - ep.2d.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "ep-2d",
                                             time = total.time["elapsed"])
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "ep-2d",
                                             lppd = lppd(X.1.test, X.2.test, y.test, ep.2d.samples))
} else if (method == "gvb") {
  library(rstan)
  mcmc.rstan <- stan_model("Hetero/Methods/MCMC.stan")
  
  start.time <- proc.time()
  
  opath <- opt_path_stan_parallel(seed_init = (seed - 1)*length(num.cores) + 1:num.cores, 
                                  seed_list = (seed - 1)*length(num.cores) + 1:num.cores, 
                                  mc.cores = num.cores, 
                                  model = mcmc.rstan,
                                  data = list(N = n,
                                              p_1 = p.1,
                                              p_2 = p.2,
                                              X_1 = X.1,
                                              X_2 = X.2,
                                              y = y,
                                              mu_theta = mu.theta,
                                              Sigma_theta = Sigma.theta),
                                  N_sam = round(mcmc.s.iter/num.cores))
  
  gvb.samples <- t(Imp_Resam_WR(opath, n_sam = mcmc.s.iter, seed = seed))
  gvb.mu <- colMeans(gvb.samples)
  gvb.Sigma <- var(gvb.samples)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = 4,
                                           method = "gvb",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], 
                                                          abs(mcmc.g.values[j, ] - dnorm(grid.points[j, ], gvb.mu[j], sqrt(gvb.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = 4,
                                                                 method = "gvb",
                                                                 mmd = max(kmmd(tail(gvb.samples, eval.size), 
                                                                                tail.mcmc.g.samples)@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = 4,
                                                     method = "gvb",
                                                     cov_norm = norm(mcmc.g.Sigma - gvb.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "gvb",
                                             time = total.time["elapsed"])
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = 4,
                                             method = "gvb",
                                             lppd = lppd(X.1.test, X.2.test, y.test, tail(gvb.samples, eval.size)))
} else {
  stop("method must be in one of mcmc, mcmc-s, ep, ep-2d, or gvb")
}

save(bench.l1.df, bench.mmd.df, bench.lppd.df, bench.cov.norm.df, bench.time.df,
     file = paste0("Hetero/Results/Big-results-", toupper(method), "-", str_pad(seed, 2, pad = "0"), ".RData"))
