
# Benchmarks for lasso linear regression

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
set.seed(seed)

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

bench.r.hat.df <- data.frame(seed = integer(),
                             bench = integer(),
                             method = character(),
                             j = integer(),
                             r_hat = double())

for (type.iter in 1:num.bench) {
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Lasso/Data/Benchmarks/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  train.ind <- sample(1:n)[1:ceiling(train.size*n)]
  X.train <- X[train.ind, , drop = F]
  y.train <- y[train.ind]
  n.train <- nrow(X.train)
  
  X.test <- X[-train.ind, , drop = F]
  y.test <- y[-train.ind]
  n.test <- nrow(X.test)
  
  ### MCMC-G
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.iter,
                   warmup = mcmc.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.g.samples <- rstan::extract(stan.res)$theta
  mcmc.g.mu <- colMeans(mcmc.g.samples)
  mcmc.g.Sigma <- var(mcmc.g.samples)
  mcmc.g.summary <- summary(stan.res)$summary
  
  grid.points <- matrix(nrow = p + 1, ncol = total.grid.points)
  mcmc.g.values <- matrix(nrow = p + 1, ncol = total.grid.points)
  
  for (j in 1:(p + 1)) {
    density.res <- density(mcmc.g.samples[, j], bw = "SJ-ste",
                           from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           n = total.grid.points)
    
    grid.points[j, ] <- density.res$x
    mcmc.g.values[j, ] <- density.res$y
  }
  
  ### MCMC
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.iter,
                   warmup = mcmc.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.samples <- rstan::extract(stan.res)$theta
  mcmc.mu <- colMeans(mcmc.samples)
  mcmc.Sigma <- var(mcmc.samples)
  mcmc.summary <- summary(stan.res)$summary
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    density.res <- density(mcmc.samples[, j], bw = "SJ-ste",
                           from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           n = total.grid.points)
    
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = type.iter,
                                           method = "mcmc",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], abs(mcmc.g.values[j, ] - density.res$y))/2)
    
    bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc",
                                                 j = j,
                                                 r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = type.iter,
                                                                 method = "mcmc",
                                                                 mmd = max(kmmd(tail(mcmc.samples, min(mcmc.iter - mcmc.warmup, eval.size)), 
                                                                                tail(mcmc.g.samples, eval.size))@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = type.iter,
                                                     method = "mcmc",
                                                     cov_norm = norm(mcmc.g.Sigma - mcmc.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "mcmc",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n.train,
                               p = p,
                               X = X.train,
                               y = y.train,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.iter,
                   warmup = mcmc.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.samples <- rstan::extract(stan.res)$theta
  
  bench.lppd.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc",
                                                 lppd = lppd(X.test, y.test, tail(mcmc.samples, min(mcmc.iter - mcmc.warmup, eval.size))))
  
  ### MCMC-A
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.a.iter,
                   warmup = mcmc.a.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.a.samples <- rstan::extract(stan.res)$theta
  mcmc.a.mu <- colMeans(mcmc.a.samples)
  mcmc.a.Sigma <- var(mcmc.a.samples)
  mcmc.a.summary <- summary(stan.res)$summary
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    density.res <- density(mcmc.a.samples[, j], bw = "SJ-ste",
                           from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           n = total.grid.points)
    
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = type.iter,
                                           method = "mcmc-a",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], abs(mcmc.g.values[j, ] - density.res$y))/2)
    
    bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc-a",
                                                 j = j,
                                                 r_hat = mcmc.a.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = type.iter,
                                                                 method = "mcmc-a",
                                                                 mmd = max(kmmd(tail(mcmc.a.samples, min(mcmc.a.iter - mcmc.a.warmup, eval.size)), 
                                                                                tail(mcmc.g.samples, eval.size))@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = type.iter,
                                                     method = "mcmc-a",
                                                     cov_norm = norm(mcmc.g.Sigma - mcmc.a.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "mcmc-a",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n.train,
                               p = p,
                               X = X.train,
                               y = y.train,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.a.iter,
                   warmup = mcmc.a.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.a.samples <- rstan::extract(stan.res)$theta
  
  bench.lppd.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc-a",
                                                 lppd = lppd(X.test, y.test, tail(mcmc.a.samples, min(mcmc.a.iter - mcmc.a.warmup, eval.size))))
  
  ### MCMC-B
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.b.iter,
                   warmup = mcmc.b.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.b.samples <- rstan::extract(stan.res)$theta
  mcmc.b.mu <- colMeans(mcmc.b.samples)
  mcmc.b.Sigma <- var(mcmc.b.samples)
  mcmc.b.summary <- summary(stan.res)$summary
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    density.res <- density(mcmc.b.samples[, j], bw = "SJ-ste",
                           from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           n = total.grid.points)
    
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = type.iter,
                                           method = "mcmc-b",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], abs(mcmc.g.values[j, ] - density.res$y))/2)
    
    bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc-b",
                                                 j = j,
                                                 r_hat = mcmc.b.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = type.iter,
                                                                 method = "mcmc-b",
                                                                 mmd = max(kmmd(tail(mcmc.b.samples, min(mcmc.b.iter - mcmc.b.warmup, eval.size)), 
                                                                                tail(mcmc.g.samples, eval.size))@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = type.iter,
                                                     method = "mcmc-b",
                                                     cov_norm = norm(mcmc.g.Sigma - mcmc.b.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "mcmc-b",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n.train,
                               p = p,
                               X = X.train,
                               y = y.train,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.b.iter,
                   warmup = mcmc.b.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.b.samples <- rstan::extract(stan.res)$theta
  
  bench.lppd.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc-b",
                                                 lppd = lppd(X.test, y.test, tail(mcmc.b.samples, min(mcmc.b.iter - mcmc.b.warmup, eval.size))))
  
  ### MCMC-C
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.c.iter,
                   warmup = mcmc.c.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.c.samples <- rstan::extract(stan.res)$theta
  mcmc.c.mu <- colMeans(mcmc.c.samples)
  mcmc.c.Sigma <- var(mcmc.c.samples)
  mcmc.c.summary <- summary(stan.res)$summary
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    density.res <- density(mcmc.c.samples[, j], bw = "SJ-ste",
                           from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                           n = total.grid.points)
    
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = type.iter,
                                           method = "mcmc-c",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], abs(mcmc.g.values[j, ] - density.res$y))/2)
    
    bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc-c",
                                                 j = j,
                                                 r_hat = mcmc.c.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = type.iter,
                                                                 method = "mcmc-c",
                                                                 mmd = max(kmmd(tail(mcmc.c.samples, min(mcmc.c.iter - mcmc.c.warmup, eval.size)), 
                                                                                tail(mcmc.g.samples, eval.size))@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = type.iter,
                                                     method = "mcmc-c",
                                                     cov_norm = norm(mcmc.g.Sigma - mcmc.c.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "mcmc-c",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                   data = list(N = n.train,
                               p = p,
                               X = X.train,
                               y = y.train,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = num.cores,
                   cores = num.cores,
                   iter = mcmc.c.iter,
                   warmup = mcmc.c.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.c.samples <- rstan::extract(stan.res)$theta
  
  bench.lppd.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 method = "mcmc-c",
                                                 lppd = lppd(X.test, y.test, tail(mcmc.c.samples, min(mcmc.c.iter - mcmc.c.warmup, eval.size))))
  
  ### EP
  
  start.time <- proc.time()
  
  ep.res <- ep(X, y, sigma.2.kappa, mu.kappa,
               lambda, eta = 0.5, alpha = 0.8, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2),
               min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = type.iter,
                                           method = "ep",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], 
                                                          abs(mcmc.g.values[j, ] - dnorm(grid.points[j, ], ep.mu[j], sqrt(ep.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = type.iter,
                                                                 method = "ep",
                                                                 mmd = max(kmmd(ep.samples, 
                                                                                tail(mcmc.g.samples, eval.size))@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = type.iter,
                                                     method = "ep",
                                                     cov_norm = norm(mcmc.g.Sigma - ep.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "ep",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ep.res <- ep(X.train, y.train, sigma.2.kappa, mu.kappa,
               lambda, eta = 0.5, alpha = 0.8, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2),
               min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "ep",
                                             fold = fold,
                                             lppd = lppd(X.test, y.test, ep.samples))
  
  ### EP-2D
  
  start.time <- proc.time()
  
  ep.2d.res <- ep_2d(X, y, sigma.2.kappa, mu.kappa,
                     lambda, eta = 0.5, alpha = 0.8, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2),
                     min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
  ep.2d.mu <- ep.2d.res$mu
  ep.2d.Sigma <- ep.2d.res$Sigma
  ep.2d.samples <- rmvnorm(eval.size, ep.2d.mu, ep.2d.Sigma)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = type.iter,
                                           method = "ep-2d",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], 
                                                          abs(mcmc.g.values[j, ] - dnorm(grid.points[j, ], ep.2d.mu[j], sqrt(ep.2d.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = type.iter,
                                                                 method = "ep-2d",
                                                                 mmd = max(kmmd(ep.2d.samples, 
                                                                                tail(mcmc.g.samples, eval.size))@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = type.iter,
                                                     method = "ep-2d",
                                                     cov_norm = norm(mcmc.g.Sigma - ep.2d.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "ep-2d",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ep.2d.res <- ep_2d(X.train, y.train, sigma.2.kappa, mu.kappa,
                     lambda, eta = 0.5, alpha = 0.8, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2),
                     min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
  ep.2d.mu <- ep.2d.res$mu
  ep.2d.Sigma <- ep.2d.res$Sigma
  ep.2d.samples <- rmvnorm(eval.size, ep.2d.mu, ep.2d.Sigma)
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "ep-2d",
                                             fold = fold,
                                             lppd = lppd(X.test, y.test, ep.2d.samples))
  
  ### MFVB
  
  start.time <- proc.time()
  
  mfvb.res <- mfvb(X, y, sigma.2.kappa, mu.kappa, lambda, maxit = 2000, tol = 1.0E-10)
  mfvb.mu <- mfvb.res$mu
  mfvb.Sigma <- mfvb.res$Sigma
  mfvb.samples <- rmvnorm(eval.size, mfvb.mu, mfvb.Sigma)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                           bench = type.iter,
                                           method = "mfvb",
                                           j = j,
                                           l1 = 1 - trapz(grid.points[j, ], 
                                                          abs(mcmc.g.values[j, ] - dnorm(grid.points[j, ], mfvb.mu[j], sqrt(mfvb.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                 bench = type.iter,
                                                                 method = "mfvb",
                                                                 mmd = max(kmmd(mfvb.samples, 
                                                                                tail(mcmc.g.samples, eval.size))@mmdstats[2], 0)))
  
  bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                     bench = type.iter,
                                                     method = "mfvb",
                                                     cov_norm = norm(mcmc.g.Sigma - mfvb.Sigma, "F"))
  
  bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "mfvb",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  mfvb.res <- mfvb(X.train, y.train, sigma.2.kappa, mu.kappa, lambda, maxit = 2000, tol = 1.0E-10)
  mfvb.mu <- mfvb.res$mu
  mfvb.Sigma <- mfvb.res$Sigma
  mfvb.samples <- rmvnorm(eval.size, mfvb.mu, mfvb.Sigma)
  
  bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "mfvb",
                                             fold = fold,
                                             lppd = lppd(X.test, y.test, mfvb.samples))
}

save(bench.l1.df, bench.mmd.df, bench.lppd.df, bench.cov.norm.df, bench.time.df, bench.r.hat.df,
     file = "Lasso/Results/Benchmarks-results.RData")
