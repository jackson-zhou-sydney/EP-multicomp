
# Simulations/benchmarks for heteroscedastic linear regression

source("EP-general-auxiliaries.R")
source("Heteroscedastic/Heteroscedastic-auxiliaries.R")

set.seed(1)

## Simulations

sim.res.df.1 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           mmd = double())

sim.res.df.2 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           fold = double(),
                           lppd = double())

sim.res.df.3 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           time = double())

sim.res.df.4 <- data.frame(sim = integer(),
                           iteration = integer(),
                           j = double(),
                           r_hat = double())

for (type.iter in 1:num.each.type) {
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    load(paste0("Heteroscedastic/Heteroscedastic-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X.1)
    p.1 <- ncol(X.1)
    p.2 <- ncol(X.2)
    mu.theta <- rep(0, p.1 + p.2)
    Sigma.theta <- sigma.2.theta*diag(p.1 + p.2)
    ind <- sample(rep(1:n.folds, ceiling(n/n.folds))[1:n])
    
    ### MCMC
    
    start.time <- proc.time()
    
    stan.res <- stan(file = "Heteroscedastic/Heteroscedastic-model.stan",
                     data = list(N = n,
                                 p_1 = p.1,
                                 p_2 = p.2,
                                 X_1 = X.1,
                                 X_2 = X.2,
                                 y = y,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta),
                     chains = mcmc.chains,
                     iter = mcmc.iter,
                     warmup = mcmc.warmup,
                     refresh = 0,
                     init = rep(0, p.1 + p.2))
    
    mcmc.samples <- rstan::extract(stan.res)$theta
    mcmc.summary <- summary(stan.res)$summary
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p.1 + p.2)) {
      sim.res.df.4 <- sim.res.df.4 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               j = j,
                                               r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
    }
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mcmc",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### MCMC-short
    
    start.time <- proc.time()
    
    stan.res <- stan(file = "Heteroscedastic/Heteroscedastic-model.stan",
                     data = list(N = n,
                                 p_1 = p.1,
                                 p_2 = p.2,
                                 X_1 = X.1,
                                 X_2 = X.2,
                                 y = y,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta),
                     chains = mcmc.chains,
                     iter = mcmc.short.iter,
                     warmup = mcmc.short.warmup,
                     refresh = 0,
                     init = rep(0, p.1 + p.2))
    
    mcmc.short.samples <- rstan::extract(stan.res)$theta
    
    total.time <- proc.time() - start.time
    
    out <- capture.output(sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "mcmc-short",
                                                                   mmd = max(kmmd(tail(mcmc.short.samples, 80), 
                                                                                  tail(mcmc.samples, 80))@mmdstats[2], 0)))
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mcmc-short",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### EP
    
    start.time <- proc.time()
    
    ep.res <- ep_c(X.1, X.2, y, Sigma.theta, mu.theta,
                   eta = 0.5, alpha = 0.75, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p.1 + p.2, p.1 + p.2),
                   min_passes = 6, max_passes = 200, tol = Inf, stop = Inf , 
                   abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    ep.samples <- rmvnorm(80, ep.mu, ep.Sigma)
    
    total.time <- proc.time() - start.time
    
    out <- capture.output(sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "ep",
                                                                   mmd = max(kmmd(ep.samples, 
                                                                                  tail(mcmc.samples, 80))@mmdstats[2], 0)))
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "ep",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### Laplace
    
    start.time <- proc.time()
    
    laplace.res <- laplace_c(X.1, X.2, y, Sigma.theta, mu.theta, rep(0, p.1 + p.2), 20000)
    laplace.mu <- laplace.res$mu
    laplace.Sigma <- laplace.res$Sigma
    laplace.samples <- rmvnorm(80, laplace.mu, laplace.Sigma)
    
    total.time <- proc.time() - start.time
    
    out <- capture.output(sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                                                   iteration = iteration,
                                                                   method = "laplace",
                                                                   mmd = max(kmmd(laplace.samples, 
                                                                                  tail(mcmc.samples, 80))@mmdstats[2], 0)))
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "laplace",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### Evaluating predictive accuracy
    
    for (fold in 1:n.folds) {
      X.1.train <- X.1[ind != fold, ]
      X.2.train <- X.2[ind != fold, ]
      y.train <- y[ind != fold]
      n.train <- nrow(X.1.train)
      
      X.1.test <- X.1[ind == fold, ]
      X.2.test <- X.2[ind == fold, ]
      y.test <- y[ind == fold]
      n.test <- nrow(X.1.test)
      
      #### MCMC-short
      
      stan.res <- stan(file = "Heteroscedastic/Heteroscedastic-model.stan",
                       data = list(N = n.train,
                                   p_1 = p.1,
                                   p_2 = p.2,
                                   X_1 = X.1.train,
                                   X_2 = X.2.train,
                                   y = y.train,
                                   mu_theta = mu.theta,
                                   Sigma_theta = Sigma.theta),
                       chains = mcmc.chains,
                       iter = mcmc.short.iter,
                       warmup = mcmc.short.warmup,
                       refresh = 0,
                       init = rep(0, p.1 + p.2))
      
      mcmc.short.samples <- rstan::extract(stan.res)$theta
      
      sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-short",
                                               fold = fold,
                                               lppd = lppd(X.1.test, X.2.test, y.test, tail(mcmc.short.samples, 80)))
      
      #### EP
      
      ep.res <- ep_c(X.1.train, X.2.train, y.train, Sigma.theta, mu.theta,
                     eta = 0.5, alpha = 0.75, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p.1 + p.2, p.1 + p.2),
                     min_passes = 6, max_passes = 200, tol = Inf, stop = Inf , 
                     abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
      ep.mu <- ep.res$mu
      ep.Sigma <- ep.res$Sigma
      ep.samples <- rmvnorm(80, ep.mu, ep.Sigma)
      
      sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "ep",
                                               fold = fold,
                                               lppd = lppd(X.1.test, X.2.test, y.test, ep.samples))
      
      #### Laplace
      
      laplace.res <- laplace_c(X.1.train, X.2.train, y.train, Sigma.theta, mu.theta, rep(0, p.1 + p.2), 20000)
      laplace.mu <- laplace.res$mu
      laplace.Sigma <- laplace.res$Sigma
      laplace.samples <- rmvnorm(80, laplace.mu, laplace.Sigma)
      
      sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "laplace",
                                               fold = fold,
                                               lppd = lppd(X.1.test, X.2.test, y.test, laplace.samples))
    }
  }
}

## Benchmarks

bench.res.df.1 <- data.frame(bench = integer(),
                             method = character(),
                             mmd = double())

bench.res.df.2 <- data.frame(bench = integer(),
                             method = character(),
                             fold = double(),
                             lppd = double())

bench.res.df.3 <- data.frame(bench = integer(),
                             method = character(),
                             time = double())

bench.res.df.4 <- data.frame(bench = integer(),
                             j = double(),
                             r_hat = double())

for (type.iter in 1:num.each.type) {
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Heteroscedastic/Heteroscedastic-data/Bench-", type.iter, ".RData"))
  n <- nrow(X.1)
  p.1 <- ncol(X.1)
  p.2 <- ncol(X.2)
  mu.theta <- rep(0, p.1 + p.2)
  Sigma.theta <- sigma.2.theta*diag(p.1 + p.2)
  ind <- sample(rep(1:n.folds, ceiling(n/n.folds))[1:n])
  
  ### MCMC
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Heteroscedastic/Heteroscedastic-model.stan",
                   data = list(N = n,
                               p_1 = p.1,
                               p_2 = p.2,
                               X_1 = X.1,
                               X_2 = X.2,
                               y = y,
                               mu_theta = mu.theta,
                               Sigma_theta = Sigma.theta),
                   chains = mcmc.chains,
                   iter = mcmc.iter,
                   warmup = mcmc.warmup,
                   refresh = 0,
                   init = rep(0, p.1 + p.2))
  
  mcmc.samples <- rstan::extract(stan.res)$theta
  mcmc.summary <- summary(stan.res)$summary
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.res.df.4 <- bench.res.df.4 %>% add_row(bench = type.iter,
                                                 j = j,
                                                 r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
  }
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "mcmc",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### MCMC-short
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Heteroscedastic/Heteroscedastic-model.stan",
                   data = list(N = n,
                               p_1 = p.1,
                               p_2 = p.2,
                               X_1 = X.1,
                               X_2 = X.2,
                               y = y,
                               mu_theta = mu.theta,
                               Sigma_theta = Sigma.theta),
                   chains = mcmc.chains,
                   iter = mcmc.short.iter,
                   warmup = mcmc.short.warmup,
                   refresh = 0,
                   init = rep(0, p.1 + p.2))
  
  mcmc.short.samples <- rstan::extract(stan.res)$theta
  
  total.time <- proc.time() - start.time
  
  out <- capture.output(bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                                     method = "mcmc-short",
                                                                     mmd = max(kmmd(tail(mcmc.short.samples, 80), 
                                                                                    tail(mcmc.samples, 80))@mmdstats[2], 0)))
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "mcmc-short",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### EP
  
  start.time <- proc.time()
  
  ep.res <- ep_c(X.1, X.2, y, Sigma.theta, mu.theta,
                 eta = 0.5, alpha = 0.75, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p.1 + p.2, p.1 + p.2),
                 min_passes = 6, max_passes = 200, tol = Inf, stop = Inf , 
                 abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  ep.samples <- rmvnorm(80, ep.mu, ep.Sigma)
  
  total.time <- proc.time() - start.time
  
  out <- capture.output(bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                                     method = "ep",
                                                                     mmd = max(kmmd(ep.samples, 
                                                                                    tail(mcmc.samples, 80))@mmdstats[2], 0)))
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "ep",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### Laplace
  
  start.time <- proc.time()
  
  laplace.res <- laplace_c(X.1, X.2, y, Sigma.theta, mu.theta, rep(0, p.1 + p.2), 20000)
  laplace.mu <- laplace.res$mu
  laplace.Sigma <- laplace.res$Sigma
  laplace.samples <- rmvnorm(80, laplace.mu, laplace.Sigma)
  
  total.time <- proc.time() - start.time
  
  out <- capture.output(bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                                     method = "laplace",
                                                                     mmd = max(kmmd(laplace.samples, 
                                                                                    tail(mcmc.samples, 80))@mmdstats[2], 0)))
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "laplace",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### Evaluating predictive accuracy
  
  for (fold in 1:n.folds) {
    X.1.train <- X.1[ind != fold, ]
    X.2.train <- X.2[ind != fold, ]
    y.train <- y[ind != fold]
    n.train <- nrow(X.1.train)
    
    X.1.test <- X.1[ind == fold, ]
    X.2.test <- X.2[ind == fold, ]
    y.test <- y[ind == fold]
    n.test <- nrow(X.1.test)
    
    #### MCMC-short
    
    stan.res <- stan(file = "Heteroscedastic/Heteroscedastic-model.stan",
                     data = list(N = n.train,
                                 p_1 = p.1,
                                 p_2 = p.2,
                                 X_1 = X.1.train,
                                 X_2 = X.2.train,
                                 y = y.train,
                                 mu_theta = mu.theta,
                                 Sigma_theta = Sigma.theta),
                     chains = mcmc.chains,
                     iter = mcmc.short.iter,
                     warmup = mcmc.short.warmup,
                     refresh = 0,
                     init = rep(0, p.1 + p.2))
    
    mcmc.short.samples <- rstan::extract(stan.res)$theta
    
    bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                 method = "mcmc-short",
                                                 fold = fold,
                                                 lppd = lppd(X.1.test, X.2.test, y.test, tail(mcmc.short.samples, 80)))
    
    #### EP
    
    ep.res <- ep_c(X.1.train, X.2.train, y.train, Sigma.theta, mu.theta,
                   eta = 0.5, alpha = 0.75, Q_star_init = 0.01*diag(2), r_star_init = rep(0, 2), offset = matrix(0, p.1 + p.2, p.1 + p.2),
                   min_passes = 6, max_passes = 200, tol = Inf, stop = Inf , 
                   abs_thresh = 0.1, rel_thresh = 0.9, delta_limit = Inf, patience = 40)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    ep.samples <- rmvnorm(80, ep.mu, ep.Sigma)
    
    bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                 method = "ep",
                                                 fold = fold,
                                                 lppd = lppd(X.1.test, X.2.test, y.test, ep.samples))
    
    #### Laplace
    
    laplace.res <- laplace_c(X.1.train, X.2.train, y.train, Sigma.theta, mu.theta, rep(0, p.1 + p.2), 20000)
    laplace.mu <- laplace.res$mu
    laplace.Sigma <- laplace.res$Sigma
    laplace.samples <- rmvnorm(80, laplace.mu, laplace.Sigma)
    
    bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                 method = "laplace",
                                                 fold = fold,
                                                 lppd = lppd(X.1.test, X.2.test, y.test, laplace.samples))
  }
}

save(sim.res.df.1, sim.res.df.2, sim.res.df.3, sim.res.df.4,
     bench.res.df.1, bench.res.df.2, bench.res.df.3, bench.res.df.4,
     file = "Heteroscedastic/Heteroscedastic-results.RData")
