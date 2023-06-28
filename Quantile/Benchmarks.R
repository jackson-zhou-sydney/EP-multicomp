
# Benchmarks for quantile regression

source("General-auxiliaries.R")
source("Quantile/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
method <- args[1]
seed <- as.numeric(args[2])
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Quantile/Methods/MCMC.stan", stanc_options = list("O1"))

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
  print(paste0("Current benchmark: ", type.iter))
  
  load(paste0("Quantile/Data/Benchmarks/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  mu.theta <- rep(0, p + 1)
  Sigma.theta <- diag(c(rep(sigma.2.beta, p), sigma.2.kappa))
  mu.beta <- mu.theta[1:p]
  Sigma.beta <- Sigma.theta[1:p, 1:p]
  mu.kappa <- mu.theta[p + 1]
  sigma.2.kappa <- Sigma.theta[p + 1, p + 1]
  
  train.ind <- sample(1:n)[1:ceiling(train.size*n)]
  X.train <- X[train.ind, , drop = F]
  y.train <- y[train.ind]
  n.train <- nrow(X.train)
  
  X.test <- X[-train.ind, , drop = F]
  y.test <- y[-train.ind]
  n.test <- nrow(X.test)
  
  if (method == "mcmc-g") {
    stan.res <- mcmc$sample(data = list(N = n,
                                        p = p,
                                        X = X,
                                        y = y,
                                        mu_theta = mu.theta,
                                        Sigma_theta = Sigma.theta,
                                        tau = tau), 
                            seed = seed, 
                            chains = num.cores, 
                            parallel_chains = num.cores,
                            iter_sampling = mcmc.g.iter,
                            iter_warmup = mcmc.g.warmup,
                            max_treedepth = max.tree.depth,
                            init = 0,
                            refresh = 1)
    
    mcmc.g.samples <- as.matrix(stan.res$draws(format = "df"))[, 2:(1 + p + 1)]
    mcmc.g.mu <- colMeans(mcmc.g.samples)
    mcmc.g.Sigma <- var(mcmc.g.samples)
    mcmc.g.summary <- stan.res$summary()
    
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
    
    tail.mcmc.g.samples <- tail(mcmc.g.samples, eval.size)
    save(tail.mcmc.g.samples, mcmc.g.mu, mcmc.g.Sigma, grid.points, mcmc.g.values, 
         file = paste0("Quantile/Results/Benchmarks-results-MCMC-G-", type.iter, "-", str_pad(seed, 2, pad = "0"), ".RData"))
  } else if (method == "mcmc") {
    load(paste0("Quantile/Results/Benchmarks-results-MCMC-G-", type.iter, "-", str_pad(seed, 2, pad = "0"), ".RData"))
    
    start.time <- proc.time()
    
    stan.res <- mcmc$sample(data = list(N = n,
                                        p = p,
                                        X = X,
                                        y = y,
                                        mu_theta = mu.theta,
                                        Sigma_theta = Sigma.theta,
                                        tau = tau), 
                            seed = seed + 1, 
                            chains = num.cores, 
                            parallel_chains = num.cores,
                            iter_sampling = mcmc.g.iter,
                            iter_warmup = mcmc.g.warmup,
                            max_treedepth = max.tree.depth,
                            init = 0,
                            refresh = 1)
    
    mcmc.samples <- as.matrix(stan.res$draws(format = "df"))[, 2:(1 + p + 1)]
    mcmc.mu <- colMeans(mcmc.samples)
    mcmc.Sigma <- var(mcmc.samples)
    mcmc.summary <- stan.res$summary()
    
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
                                                   r_hat = mcmc.summary %>% filter(variable == paste0("theta[", j, "]")) %>% pull(rhat))
    }
    
    out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                   bench = type.iter,
                                                                   method = "mcmc",
                                                                   mmd = max(kmmd(tail(mcmc.samples, eval.size), 
                                                                                  tail.mcmc.g.samples)@mmdstats[2], 0)))
    
    bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                       bench = type.iter,
                                                       method = "mcmc",
                                                       cov_norm = norm(mcmc.g.Sigma - mcmc.Sigma, "F"))
    
    bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "mcmc",
                                               time = total.time["elapsed"])
    
    stan.res <- mcmc$sample(data = list(N = n.train,
                                        p = p,
                                        X = X.train,
                                        y = y.train,
                                        mu_theta = mu.theta,
                                        Sigma_theta = Sigma.theta,
                                        tau = tau), 
                            seed = seed + 1, 
                            chains = num.cores, 
                            parallel_chains = num.cores,
                            iter_sampling = mcmc.g.iter,
                            iter_warmup = mcmc.g.warmup,
                            max_treedepth = max.tree.depth,
                            init = 0,
                            refresh = 1)
    
    mcmc.samples <- as.matrix(stan.res$draws(format = "df"))[, 2:(1 + p + 1)]
    
    bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "mcmc",
                                               lppd = lppd(X.test, y.test, tau, tail(mcmc.samples, eval.size)))
  } else if (method == "mcmc-s") {
    load(paste0("Quantile/Results/Benchmarks-results-MCMC-G-", type.iter, "-", str_pad(seed, 2, pad = "0"), ".RData"))
    
    load("Quantile/Results/Benchmarks-conv-table.RData")
    mcmc.test.iter <- bench.r.hat.table %>% pull(mcmc_iter) %>% unique() %>% sort()
    ind <- which(mcmc.test.iter == bench.r.hat.table %>% filter(sim == type.iter) %>% filter(mean_max_r_hat > r.hat.tol) %>% pull(mcmc_iter) %>% max()) + 1
    mcmc.s.iter <- mcmc.test.iter[min(ind, length(mcmc.test.iter))]
    
    start.time <- proc.time()
    
    stan.res <- mcmc$sample(data = list(N = n,
                                        p = p,
                                        X = X,
                                        y = y,
                                        mu_theta = mu.theta,
                                        Sigma_theta = Sigma.theta,
                                        tau = tau), 
                            seed = seed + 2, 
                            chains = num.cores, 
                            parallel_chains = num.cores,
                            iter_sampling = mcmc.s.iter,
                            iter_warmup = mcmc.s.warmup,
                            max_treedepth = max.tree.depth,
                            init = 0,
                            refresh = 1)
    
    mcmc.s.samples <- as.matrix(stan.res$draws(format = "df"))[, 2:(1 + p + 1)]
    mcmc.s.mu <- colMeans(mcmc.s.samples)
    mcmc.s.Sigma <- var(mcmc.s.samples)
    mcmc.s.summary <- stan.res$summary()
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p + 1)) {
      density.res <- density(mcmc.s.samples[, j], bw = "SJ-ste",
                             from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                             to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                             n = total.grid.points)
      
      bench.l1.df <- bench.l1.df %>% add_row(seed = seed,
                                             bench = type.iter,
                                             method = "mcmc-s",
                                             j = j,
                                             l1 = 1 - trapz(grid.points[j, ], abs(mcmc.g.values[j, ] - density.res$y))/2)
      
      bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                   bench = type.iter,
                                                   method = "mcmc-s",
                                                   j = j,
                                                   r_hat = mcmc.s.summary %>% filter(variable == paste0("theta[", j, "]")) %>% pull(rhat))
    }
    
    out <- capture.output(bench.mmd.df <- bench.mmd.df %>% add_row(seed = seed,
                                                                   bench = type.iter,
                                                                   method = "mcmc-s",
                                                                   mmd = max(kmmd(tail(mcmc.s.samples, eval.size), 
                                                                                  tail.mcmc.g.samples)@mmdstats[2], 0)))
    
    bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                       bench = type.iter,
                                                       method = "mcmc-s",
                                                       cov_norm = norm(mcmc.g.Sigma - mcmc.s.Sigma, "F"))
    
    bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "mcmc-s",
                                               time = total.time["elapsed"])
    
    stan.res <- mcmc$sample(data = list(N = n.train,
                                        p = p,
                                        X = X.train,
                                        y = y.train,
                                        mu_theta = mu.theta,
                                        Sigma_theta = Sigma.theta,
                                        tau = tau), 
                            seed = seed + 2, 
                            chains = num.cores, 
                            parallel_chains = num.cores,
                            iter_sampling = mcmc.s.iter,
                            iter_warmup = mcmc.s.warmup,
                            max_treedepth = max.tree.depth,
                            init = 0,
                            refresh = 1)
    
    mcmc.s.samples <- as.matrix(stan.res$draws(format = "df"))[, 2:(1 + p + 1)]
    
    bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "mcmc-s",
                                               lppd = lppd(X.test, y.test, tau, tail(mcmc.s.samples, eval.size)))
  } else if (method == "ep") {
    load(paste0("Quantile/Results/Benchmarks-results-MCMC-G-", type.iter, "-", str_pad(seed, 2, pad = "0"), ".RData"))
    
    start.time <- proc.time()
    
    ep.res <- ep(X, y, Sigma.theta, mu.theta,
                 tau, eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
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
                                                                                  tail.mcmc.g.samples)@mmdstats[2], 0)))
    
    bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                       bench = type.iter,
                                                       method = "ep",
                                                       cov_norm = norm(mcmc.g.Sigma - ep.Sigma, "F"))
    
    bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "ep",
                                               time = total.time["elapsed"])
    
    ep.res <- ep(X.train, y.train, Sigma.theta, mu.theta,
                 tau, eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
                 min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
    
    bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "ep",
                                               lppd = lppd(X.test, y.test, tau, ep.samples))
  } else if (method == "ep-2d") {
    load(paste0("Quantile/Results/Benchmarks-results-MCMC-G-", type.iter, "-", str_pad(seed, 2, pad = "0"), ".RData"))
    
    start.time <- proc.time()
    
    ep.2d.res <- ep_2d(X, y, Sigma.theta, mu.theta,
                       tau, eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
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
                                                                                  tail.mcmc.g.samples)@mmdstats[2], 0)))
    
    bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                       bench = type.iter,
                                                       method = "ep-2d",
                                                       cov_norm = norm(mcmc.g.Sigma - ep.2d.Sigma, "F"))
    
    bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "ep-2d",
                                               time = total.time["elapsed"])
    
    ep.2d.res <- ep_2d(X.train, y.train, Sigma.theta, mu.theta,
                       tau, eta = 0.5, alpha = 0.5, Q_star_init = diag(2), r_star_init = rep(0, 2),
                       min_passes = 6, max_passes = 200, thresh = 0.05, n_grid = 400, verbose = F)
    ep.2d.mu <- ep.2d.res$mu
    ep.2d.Sigma <- ep.2d.res$Sigma
    ep.2d.samples <- rmvnorm(eval.size, ep.2d.mu, ep.2d.Sigma)
    
    bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "ep-2d",
                                               lppd = lppd(X.test, y.test, tau, ep.2d.samples))
  } else if (method == "mfvb") {
    load(paste0("Quantile/Results/Benchmarks-results-MCMC-G-", type.iter, "-", str_pad(seed, 2, pad = "0"), ".RData"))
    
    start.time <- proc.time()
    
    mfvb.res <- mfvb(X, y, Sigma.beta, mu.beta, sigma.2.kappa, mu.kappa, tau,
                     min_iter = 6, max_iter = 200, thresh = 0.05, n_grid = 400, verbose = F)
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
                                                                                  tail.mcmc.g.samples)@mmdstats[2], 0)))
    
    bench.cov.norm.df <- bench.cov.norm.df %>% add_row(seed = seed,
                                                       bench = type.iter,
                                                       method = "mfvb",
                                                       cov_norm = norm(mcmc.g.Sigma - mfvb.Sigma, "F"))
    
    bench.time.df <- bench.time.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "mfvb",
                                               time = total.time["elapsed"])
    
    mfvb.res <- mfvb(X.train, y.train, Sigma.beta, mu.beta, sigma.2.kappa, mu.kappa, tau,
                     min_iter = 6, max_iter = 200, thresh = 0.05, n_grid = 400, verbose = F)
    mfvb.mu <- mfvb.res$mu
    mfvb.Sigma <- mfvb.res$Sigma
    mfvb.samples <- rmvnorm(eval.size, mfvb.mu, mfvb.Sigma)
    
    bench.lppd.df <- bench.lppd.df %>% add_row(seed = seed,
                                               bench = type.iter,
                                               method = "mfvb",
                                               lppd = lppd(X.test, y.test, tau, mfvb.samples))
  } else {
    stop("method must be in one of mcmc-g, mcmc, mcmc-s, ep, ep-2d, or mfvb")
  }
}

if (method != "mcmc-g") {
  save(bench.l1.df, bench.mmd.df, bench.lppd.df, bench.cov.norm.df, bench.time.df, bench.r.hat.df,
       file = paste0("Quantile/Results/Benchmarks-results-", toupper(method), "-", str_pad(seed, 2, pad = "0"), ".RData"))
}
