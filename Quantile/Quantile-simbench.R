
# Simulations/benchmarks for quantile regression

source("EP-general-auxiliaries.R")
source("Quantile/Quantile-auxiliaries.R")

set.seed(1)

## Simulations

sim.res.df.1 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           j = integer(),
                           l1 = double())

sim.res.df.2 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           mmd = double(),
                           cov_norm = double())

sim.res.df.3 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           fold = integer(),
                           lppd = double())

sim.res.df.4 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           time = double())

sim.res.df.5 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           j = integer(),
                           r_hat = double())

for (type.iter in 1:num.each.type) {
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    load(paste0("Quantile/Quantile-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X)
    p <- ncol(X)
    mu.theta <- rep(0, p + 1)
    Sigma.theta <- sigma.2.theta*diag(p + 1)
    ind <- sample(rep(1:n.folds, ceiling(n/n.folds))[1:n])
    
    mu.beta <- mu.theta[1:p]
    Sigma.beta <- Sigma.theta[1:p, 1:p]
    mu.kappa <- mu.theta[p + 1]
    sigma.2.kappa <- Sigma.theta[p + 1, p + 1]
    
    ### MCMC
    
    start.time <- proc.time()
    
    stan.res <- stan(file = "Quantile/Quantile-model.stan",
                     data = list(N = n,
                                 p = p,
                                 X = X,
                                 y = y,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta,
                                 tau = tau),
                     chains = mcmc.chains,
                     iter = mcmc.iter,
                     warmup = mcmc.warmup,
                     refresh = 0,
                     init = rep(0, p + 1))
    
    mcmc.samples <- rstan::extract(stan.res)$theta
    mcmc.mu <- colMeans(mcmc.samples)
    mcmc.Sigma <- var(mcmc.samples)
    mcmc.summary <- summary(stan.res)$summary
    
    total.time <- proc.time() - start.time
    
    grid.points <- matrix(nrow = p + 1, ncol = total.grid.points)
    mcmc.values <- matrix(nrow = p + 1, ncol = total.grid.points)
    
    for (j in 1:(p + 1)) {
      density.res <- density(mcmc.samples[, j], bw = "SJ-ste",
                             from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             n = total.grid.points)
      
      grid.points[j, ] <- density.res$x
      mcmc.values[j, ] <- density.res$y
      
      sim.res.df.5 <- sim.res.df.5 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc",
                                               j = j,
                                               r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
    }
    
    sim.res.df.4 <- sim.res.df.4 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mcmc",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### MCMC-A
    
    start.time <- proc.time()
    
    stan.res <- stan(file = "Quantile/Quantile-model.stan",
                     data = list(N = n,
                                 p = p,
                                 X = X,
                                 y = y,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta,
                                 tau = tau),
                     chains = mcmc.chains,
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
                             from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             n = total.grid.points)
      
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-a",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - density.res$y))/2)
      
      sim.res.df.5 <- sim.res.df.5 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-a",
                                               j = j,
                                               r_hat = mcmc.a.summary[paste0("theta[", j, "]"), "Rhat"])
    }
    
    out <- capture.output(sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "mcmc-a",
                                                                   mmd = max(kmmd(tail(mcmc.a.samples, min(mcmc.a.iter - mcmc.a.warmup, eval.size)), 
                                                                                  tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                   cov_norm = norm(mcmc.Sigma - mcmc.a.Sigma, "F")))
    
    sim.res.df.4 <- sim.res.df.4 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mcmc-a",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### MCMC-B
    
    start.time <- proc.time()
    
    stan.res <- stan(file = "Quantile/Quantile-model.stan",
                     data = list(N = n,
                                 p = p,
                                 X = X,
                                 y = y,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta,
                                 tau = tau),
                     chains = mcmc.chains,
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
                             from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             n = total.grid.points)
      
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-b",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - density.res$y))/2)
      
      sim.res.df.5 <- sim.res.df.5 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-b",
                                               j = j,
                                               r_hat = mcmc.b.summary[paste0("theta[", j, "]"), "Rhat"])
    }
    
    out <- capture.output(sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "mcmc-b",
                                                                   mmd = max(kmmd(tail(mcmc.b.samples, min(mcmc.b.iter - mcmc.b.warmup, eval.size)), 
                                                                                  tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                   cov_norm = norm(mcmc.Sigma - mcmc.b.Sigma, "F")))
    
    sim.res.df.4 <- sim.res.df.4 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mcmc-b",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### MCMC-C
    
    start.time <- proc.time()
    
    stan.res <- stan(file = "Quantile/Quantile-model.stan",
                     data = list(N = n,
                                 p = p,
                                 X = X,
                                 y = y,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta,
                                 tau = tau),
                     chains = mcmc.chains,
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
                             from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                             n = total.grid.points)
      
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-c",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - density.res$y))/2)
      
      sim.res.df.5 <- sim.res.df.5 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-c",
                                               j = j,
                                               r_hat = mcmc.c.summary[paste0("theta[", j, "]"), "Rhat"])
    }
    
    out <- capture.output(sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "mcmc-c",
                                                                   mmd = max(kmmd(tail(mcmc.c.samples, min(mcmc.c.iter - mcmc.c.warmup, eval.size)), 
                                                                                  tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                   cov_norm = norm(mcmc.Sigma - mcmc.c.Sigma, "F")))
    
    sim.res.df.4 <- sim.res.df.4 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mcmc-c",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### EP
    
    start.time <- proc.time()
    
    ep.res <- ep_c(X, y, Sigma.theta, mu.theta,
                   tau, eta = 0.5, alpha = 1, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p + 1, p + 1),
                   min_passes = 6, max_passes = 200, tol = Inf, stop = Inf,
                   abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p + 1)) {
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "ep",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], 
                                                              abs(mcmc.values[j, ] - dnorm(grid.points[j, ], ep.mu[j], sqrt(ep.Sigma[j, j]))))/2)
    }
    
    out <- capture.output(sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "ep",
                                                                   mmd = max(kmmd(ep.samples, 
                                                                                  tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                   cov_norm = norm(mcmc.Sigma - ep.Sigma, "F")))
    
    sim.res.df.4 <- sim.res.df.4 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "ep",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### MFVB
    
    start.time <- proc.time()
    
    mfvb.res <- mfvb_c(X, y, Sigma.beta, mu.beta, sigma.2.kappa, mu.kappa,
                       tau, maxit = 2000, tol = 1.0E-10)
    mfvb.mu <- mfvb.res$mu
    mfvb.Sigma <- mfvb.res$Sigma
    mfvb.samples <- rmvnorm(eval.size, mfvb.mu, mfvb.Sigma)
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p + 1)) {
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mfvb",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], 
                                                              abs(mcmc.values[j, ] - dnorm(grid.points[j, ], mfvb.mu[j], sqrt(mfvb.Sigma[j, j]))))/2)
    }
    
    out <- capture.output(sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "mfvb",
                                                                   mmd = max(kmmd(mfvb.samples, 
                                                                                  tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                   cov_norm = norm(mcmc.Sigma - mfvb.Sigma, "F")))
    
    sim.res.df.4 <- sim.res.df.4 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mfvb",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### Evaluating predictive accuracy
    
    for (fold in 1:n.folds) {
      X.train <- X[ind != fold, , drop = F]
      y.train <- y[ind != fold]
      n.train <- nrow(X.train)
      
      X.test <- X[ind == fold, , drop = F]
      y.test <- y[ind == fold]
      n.test <- nrow(X.test)
      
      #### MCMC-A
      
      stan.res <- stan(file = "Quantile/Quantile-model.stan",
                       data = list(N = n.train,
                                   p = p,
                                   X = X.train,
                                   y = y.train,
                                   mu_theta = mu.theta,
                                   Sigma_theta = Sigma.theta,
                                   tau = tau),
                       chains = mcmc.chains,
                       iter = mcmc.a.iter,
                       warmup = mcmc.a.warmup,
                       refresh = 0,
                       init = rep(0, p + 1))
      
      mcmc.a.samples <- rstan::extract(stan.res)$theta
      
      sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-a",
                                               fold = fold,
                                               lppd = lppd(X.test, y.test, tau, tail(mcmc.a.samples, min(mcmc.a.iter - mcmc.a.warmup, eval.size))))
      
      #### MCMC-B
      
      stan.res <- stan(file = "Quantile/Quantile-model.stan",
                       data = list(N = n.train,
                                   p = p,
                                   X = X.train,
                                   y = y.train,
                                   mu_theta = mu.theta,
                                   Sigma_theta = Sigma.theta,
                                   tau = tau),
                       chains = mcmc.chains,
                       iter = mcmc.b.iter,
                       warmup = mcmc.b.warmup,
                       refresh = 0,
                       init = rep(0, p + 1))
      
      mcmc.b.samples <- rstan::extract(stan.res)$theta
      
      sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-b",
                                               fold = fold,
                                               lppd = lppd(X.test, y.test, tau, tail(mcmc.b.samples, min(mcmc.b.iter - mcmc.b.warmup, eval.size))))
      
      #### MCMC-C
      
      stan.res <- stan(file = "Quantile/Quantile-model.stan",
                       data = list(N = n.train,
                                   p = p,
                                   X = X.train,
                                   y = y.train,
                                   mu_theta = mu.theta,
                                   Sigma_theta = Sigma.theta,
                                   tau = tau),
                       chains = mcmc.chains,
                       iter = mcmc.c.iter,
                       warmup = mcmc.c.warmup,
                       refresh = 0,
                       init = rep(0, p + 1))
      
      mcmc.c.samples <- rstan::extract(stan.res)$theta
      
      sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-c",
                                               fold = fold,
                                               lppd = lppd(X.test, y.test, tau, tail(mcmc.c.samples, min(mcmc.c.iter - mcmc.c.warmup, eval.size))))
      
      #### EP
      
      ep.res <- ep_c(X.train, y.train, Sigma.theta, mu.theta,
                     tau, eta = 0.5, alpha = 1, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p + 1, p + 1),
                     min_passes = 6, max_passes = 200, tol = Inf, stop = Inf,
                     abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
      ep.mu <- ep.res$mu
      ep.Sigma <- ep.res$Sigma
      ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
      
      sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "ep",
                                               fold = fold,
                                               lppd = lppd(X.test, y.test, tau, ep.samples))
      
      #### MFVB
      
      mfvb.res <- mfvb_c(X.train, y.train, Sigma.beta, mu.beta, sigma.2.kappa, mu.kappa,
                         tau, maxit = 2000, tol = 1.0E-10)
      mfvb.mu <- mfvb.res$mu
      mfvb.Sigma <- mfvb.res$Sigma
      mfvb.samples <- rmvnorm(eval.size, mfvb.mu, mfvb.Sigma)
      
      sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mfvb",
                                               fold = fold,
                                               lppd = lppd(X.test, y.test, tau, mfvb.samples))
    }
  }
}

## Benchmarks

bench.res.df.1 <- data.frame(bench = integer(),
                             method = character(),
                             j = integer(),
                             l1 = double())

bench.res.df.2 <- data.frame(bench = integer(),
                             method = character(),
                             mmd = double(),
                             cov_norm = double())

bench.res.df.3 <- data.frame(bench = integer(),
                             method = character(),
                             fold = integer(),
                             lppd = double())

bench.res.df.4 <- data.frame(bench = integer(),
                             method = character(),
                             time = double())

bench.res.df.5 <- data.frame(bench = integer(),
                             method = character(),
                             j = integer(),
                             r_hat = double())

for (type.iter in 1:num.each.type) {
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Quantile/Quantile-data/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  mu.theta <- rep(0, p + 1)
  Sigma.theta <- sigma.2.theta*diag(p + 1)
  ind <- sample(rep(1:n.folds, ceiling(n/n.folds))[1:n])
  
  mu.beta <- mu.theta[1:p]
  Sigma.beta <- Sigma.theta[1:p, 1:p]
  mu.kappa <- mu.theta[p + 1]
  sigma.2.kappa <- Sigma.theta[p + 1, p + 1]
  
  ### MCMC
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Quantile/Quantile-model.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_theta = mu.theta,
                               Sigma_theta = Sigma.theta,
                               tau = tau),
                   chains = mcmc.chains,
                   iter = mcmc.iter,
                   warmup = mcmc.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.samples <- rstan::extract(stan.res)$theta
  mcmc.mu <- colMeans(mcmc.samples)
  mcmc.Sigma <- var(mcmc.samples)
  mcmc.summary <- summary(stan.res)$summary
  
  total.time <- proc.time() - start.time
  
  grid.points <- matrix(nrow = p + 1, ncol = total.grid.points)
  mcmc.values <- matrix(nrow = p + 1, ncol = total.grid.points)
  
  for (j in 1:(p + 1)) {
    density.res <- density(mcmc.samples[, j], bw = "SJ-ste",
                           from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           n = total.grid.points)
    
    grid.points[j, ] <- density.res$x
    mcmc.values[j, ] <- density.res$y
    
    bench.res.df.5 <- bench.res.df.5 %>% add_row(bench = type.iter,
                                                 method = "mcmc",
                                                 j = j,
                                                 r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  bench.res.df.4 <- bench.res.df.4 %>% add_row(bench = type.iter,
                                               method = "mcmc",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### MCMC-A
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Quantile/Quantile-model.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_theta = mu.theta,
                               Sigma_theta = Sigma.theta,
                               tau = tau),
                   chains = mcmc.chains,
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
                           from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           n = total.grid.points)
    
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "mcmc-a",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - density.res$y))/2)
    
    bench.res.df.5 <- bench.res.df.5 %>% add_row(bench = type.iter,
                                                 method = "mcmc-a",
                                                 j = j,
                                                 r_hat = mcmc.a.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  out <- capture.output(bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                                     method = "mcmc-a",
                                                                     mmd = max(kmmd(tail(mcmc.a.samples, min(mcmc.a.iter - mcmc.a.warmup, eval.size)), 
                                                                                    tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                     cov_norm = norm(mcmc.Sigma - mcmc.a.Sigma, "F")))
  
  bench.res.df.4 <- bench.res.df.4 %>% add_row(bench = type.iter,
                                               method = "mcmc-a",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### MCMC-B
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Quantile/Quantile-model.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_theta = mu.theta,
                               Sigma_theta = Sigma.theta,
                               tau = tau),
                   chains = mcmc.chains,
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
                           from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           n = total.grid.points)
    
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "mcmc-b",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - density.res$y))/2)
    
    bench.res.df.5 <- bench.res.df.5 %>% add_row(bench = type.iter,
                                                 method = "mcmc-b",
                                                 j = j,
                                                 r_hat = mcmc.b.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  out <- capture.output(bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                                     method = "mcmc-b",
                                                                     mmd = max(kmmd(tail(mcmc.b.samples, min(mcmc.b.iter - mcmc.b.warmup, eval.size)), 
                                                                                    tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                     cov_norm = norm(mcmc.Sigma - mcmc.b.Sigma, "F")))
  
  bench.res.df.4 <- bench.res.df.4 %>% add_row(bench = type.iter,
                                               method = "mcmc-b",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### MCMC-C
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Quantile/Quantile-model.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_theta = mu.theta,
                               Sigma_theta = Sigma.theta,
                               tau = tau),
                   chains = mcmc.chains,
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
                           from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                           n = total.grid.points)
    
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "mcmc-c",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - density.res$y))/2)
    
    bench.res.df.5 <- bench.res.df.5 %>% add_row(bench = type.iter,
                                                 method = "mcmc-c",
                                                 j = j,
                                                 r_hat = mcmc.c.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  out <- capture.output(bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                                     method = "mcmc-c",
                                                                     mmd = max(kmmd(tail(mcmc.c.samples, min(mcmc.c.iter - mcmc.c.warmup, eval.size)), 
                                                                                    tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                     cov_norm = norm(mcmc.Sigma - mcmc.c.Sigma, "F")))
  
  bench.res.df.4 <- bench.res.df.4 %>% add_row(bench = type.iter,
                                               method = "mcmc-c",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### EP
  
  start.time <- proc.time()
  
  ep.res <- ep_c(X, y, Sigma.theta, mu.theta,
                 tau, eta = 0.5, alpha = 1, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p + 1, p + 1),
                 min_passes = 6, max_passes = 200, tol = Inf, stop = Inf,
                 abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "ep",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], 
                                                                abs(mcmc.values[j, ] - dnorm(grid.points[j, ], ep.mu[j], sqrt(ep.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                                     method = "ep",
                                                                     mmd = max(kmmd(ep.samples, 
                                                                                    tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                     cov_norm = norm(mcmc.Sigma - ep.Sigma, "F")))
  
  bench.res.df.4 <- bench.res.df.4 %>% add_row(bench = type.iter,
                                               method = "ep",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### MFVB
  
  start.time <- proc.time()
  
  mfvb.res <- mfvb_c(X, y, Sigma.beta, mu.beta, sigma.2.kappa, mu.kappa,
                     tau, maxit = 2000, tol = 1.0E-10)
  mfvb.mu <- mfvb.res$mu
  mfvb.Sigma <- mfvb.res$Sigma
  mfvb.samples <- rmvnorm(eval.size, mfvb.mu, mfvb.Sigma)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "mfvb",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], 
                                                                abs(mcmc.values[j, ] - dnorm(grid.points[j, ], mfvb.mu[j], sqrt(mfvb.Sigma[j, j]))))/2)
  }
  
  out <- capture.output(bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                                     method = "mfvb",
                                                                     mmd = max(kmmd(mfvb.samples, 
                                                                                    tail(mcmc.samples, eval.size))@mmdstats[2], 0),
                                                                     cov_norm = norm(mcmc.Sigma - mfvb.Sigma, "F")))
  
  bench.res.df.4 <- bench.res.df.4 %>% add_row(bench = type.iter,
                                               method = "mfvb",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### Evaluating predictive accuracy
  
  for (fold in 1:n.folds) {
    X.train <- X[ind != fold, , drop = F]
    y.train <- y[ind != fold]
    n.train <- nrow(X.train)
    
    X.test <- X[ind == fold, , drop = F]
    y.test <- y[ind == fold]
    n.test <- nrow(X.test)
    
    #### MCMC-A
    
    stan.res <- stan(file = "Quantile/Quantile-model.stan",
                     data = list(N = n.train,
                                 p = p,
                                 X = X.train,
                                 y = y.train,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta,
                                 tau = tau),
                     chains = mcmc.chains,
                     iter = mcmc.a.iter,
                     warmup = mcmc.a.warmup,
                     refresh = 0,
                     init = rep(0, p + 1))
    
    mcmc.a.samples <- rstan::extract(stan.res)$theta
    
    bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                                 method = "mcmc-a",
                                                 fold = fold,
                                                 lppd = lppd(X.test, y.test, tau, tail(mcmc.a.samples, min(mcmc.a.iter - mcmc.a.warmup, eval.size))))
    
    #### MCMC-B
    
    stan.res <- stan(file = "Quantile/Quantile-model.stan",
                     data = list(N = n.train,
                                 p = p,
                                 X = X.train,
                                 y = y.train,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta,
                                 tau = tau),
                     chains = mcmc.chains,
                     iter = mcmc.b.iter,
                     warmup = mcmc.b.warmup,
                     refresh = 0,
                     init = rep(0, p + 1))
    
    mcmc.b.samples <- rstan::extract(stan.res)$theta
    
    bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                                 method = "mcmc-b",
                                                 fold = fold,
                                                 lppd = lppd(X.test, y.test, tau, tail(mcmc.b.samples, min(mcmc.b.iter - mcmc.b.warmup, eval.size))))
    
    #### MCMC-C
    
    stan.res <- stan(file = "Quantile/Quantile-model.stan",
                     data = list(N = n.train,
                                 p = p,
                                 X = X.train,
                                 y = y.train,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta,
                                 tau = tau),
                     chains = mcmc.chains,
                     iter = mcmc.c.iter,
                     warmup = mcmc.c.warmup,
                     refresh = 0,
                     init = rep(0, p + 1))
    
    mcmc.c.samples <- rstan::extract(stan.res)$theta
    
    bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                                 method = "mcmc-c",
                                                 fold = fold,
                                                 lppd = lppd(X.test, y.test, tau, tail(mcmc.c.samples, min(mcmc.c.iter - mcmc.c.warmup, eval.size))))
    
    #### EP
    
    ep.res <- ep_c(X.train, y.train, Sigma.theta, mu.theta,
                   tau, eta = 0.5, alpha = 1, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p + 1, p + 1),
                   min_passes = 6, max_passes = 200, tol = Inf, stop = Inf,
                   abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    ep.samples <- rmvnorm(eval.size, ep.mu, ep.Sigma)
    
    bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                                 method = "ep",
                                                 fold = fold,
                                                 lppd = lppd(X.test, y.test, tau, ep.samples))
    
    #### MFVB
    
    mfvb.res <- mfvb_c(X.train, y.train, Sigma.beta, mu.beta, sigma.2.kappa, mu.kappa,
                       tau, maxit = 2000, tol = 1.0E-10)
    mfvb.mu <- mfvb.res$mu
    mfvb.Sigma <- mfvb.res$Sigma
    mfvb.samples <- rmvnorm(eval.size, mfvb.mu, mfvb.Sigma)
    
    bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                                 method = "mfvb",
                                                 fold = fold,
                                                 lppd = lppd(X.test, y.test, tau, mfvb.samples))
  }
}

save(sim.res.df.1, sim.res.df.2, sim.res.df.3, sim.res.df.4, sim.res.df.5,
     bench.res.df.1, bench.res.df.2, bench.res.df.3, bench.res.df.4, bench.res.df.5,
     file = "Quantile/Quantile-results.RData")
