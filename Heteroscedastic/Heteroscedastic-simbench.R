
# Simulations/benchmarks for heteroscedastic linear regression

source("EP-general-auxiliaries.R")
source("Heteroscedastic/Heteroscedastic-auxiliaries.R")

set.seed(1)

## Simulations

sim.res.df.1 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           j = double(),
                           l1 = double(),
                           diff_mu = double(),
                           diff_sigma = double())

sim.res.df.2 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           repetition = double(),
                           match_pairs = double())

sim.res.df.3 <- data.frame(sim = integer(),
                           iteration = integer(),
                           method = character(),
                           time = double())

sim.res.df.4 <- data.frame(sim = integer(),
                           iteration = integer(),
                           j = double(),
                           r_hat = double())

sim.res.list <- list()

for (type.iter in 1:num.each.type) {
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    load(paste0("Heteroscedastic/Heteroscedastic-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X.1)
    p.1 <- ncol(X.1)
    p.2 <- ncol(X.2)
    mu.theta <- rep(0, p.1 + p.2)
    Sigma.theta <- sigma.2.theta*diag(p.1 + p.2)
    
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
    mcmc.mu <- colMeans(mcmc.samples)
    mcmc.Sigma <- var(mcmc.samples)
    mcmc.summary <- summary(stan.res)$summary
    
    total.time <- proc.time() - start.time
    
    grid.points <- matrix(nrow = p.1 + p.2, ncol = total.grid.points)
    mcmc.values <- matrix(nrow = p.1 + p.2, ncol = total.grid.points)
    
    for (j in 1:(p.1 + p.2)) {
      grid.points[j, ] <- seq(from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                              to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                              length = total.grid.points)
      
      mcmc.values[j, ] <- demp(grid.points[j, ], obs = mcmc.samples[, j])
      
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
    mcmc.short.mu <- colMeans(mcmc.short.samples)
    mcmc.short.Sigma <- var(mcmc.short.samples)
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p.1 + p.2)) {
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mcmc-short",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - demp(grid.points[j, ], 
                                                                                                            obs = mcmc.short.samples[, j])))/2,
                                               diff_mu = (mcmc.short.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                               diff_sigma = (sqrt(mcmc.short.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
    }
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mcmc-short",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### EP
    
    start.time <- proc.time()
    
    ep.res <- ep.approx(X.1, X.2, y, mu.theta, Sigma.theta, 
                        eta = 0.5, alpha = 1, Q.star.init = 0.01*diag(2), r.star.init = rep(0, 2), offset = 0,
                        min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                        abs.thresh = 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p.1 + p.2)) {
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "ep",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - dnorm(grid.points[j, ], 
                                                                                                             mean = ep.mu[j], 
                                                                                                             sd = sqrt(ep.Sigma[j, j]))))/2,
                                               diff_mu = (ep.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                               diff_sigma = (sqrt(ep.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
    }
    
    for (repetition in 1:match.reps) {
      index <- (nrow(mcmc.samples) - match.size*repetition + 1):(nrow(mcmc.samples) - match.size*(repetition - 1))
      
      sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "ep",
                                               repetition = repetition,
                                               match_pairs = nbp.match.pairs(rmvnorm(match.size, ep.mu, ep.Sigma),
                                                                             unname(mcmc.samples[index, ])))
    }
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "ep",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### Laplace
    
    start.time <- proc.time()
    
    laplace.res <- laplace.approx(X.1, X.2, y, mu.theta, Sigma.theta, lambda.init = 0.5, maxit = 50000)
    laplace.mu <- laplace.res$mu
    laplace.Sigma <- laplace.res$Sigma
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p.1 + p.2)) {
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "laplace",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - dnorm(grid.points[j, ], 
                                                                                                             mean = laplace.mu[j], 
                                                                                                             sd = sqrt(laplace.Sigma[j, j]))))/2,
                                               diff_mu = (laplace.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                               diff_sigma = (sqrt(laplace.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
    }
    
    for (repetition in 1:match.reps) {
      index <- (nrow(mcmc.samples) - match.size*repetition + 1):(nrow(mcmc.samples) - match.size*(repetition - 1))
      
      sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "laplace",
                                               repetition = repetition,
                                               match_pairs = nbp.match.pairs(rmvnorm(match.size, laplace.mu, laplace.Sigma),
                                                                             unname(mcmc.samples[index, ])))
    }
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "laplace",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### Means and covariances
    
    sim.res.list[[paste0("s", type.iter)]][[paste0("i", iteration)]] <- list(mcmc.mu = as.vector(mcmc.mu), mcmc.Sigma = mcmc.Sigma,
                                                                             mcmc.short.mu = as.vector(mcmc.short.mu), mcmc.short.Sigma = mcmc.short.Sigma,
                                                                             ep.mu = as.vector(ep.mu), ep.Sigma = ep.Sigma,
                                                                             laplace.mu = as.vector(laplace.mu), laplace.Sigma = laplace.Sigma)
  }
}

## Benchmarks

bench.res.df.1 <- data.frame(bench = integer(),
                             method = character(),
                             j = double(),
                             l1 = double(),
                             diff_mu = double(),
                             diff_sigma = double())

bench.res.df.2 <- data.frame(bench = integer(),
                             method = character(),
                             repetition = double(),
                             match_pairs = double())

bench.res.df.3 <- data.frame(bench = integer(),
                             method = character(),
                             time = double())

bench.res.df.4 <- data.frame(bench = integer(),
                             j = double(),
                             r_hat = double())

bench.res.list <- list()

for (type.iter in 1:num.each.type) {
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Heteroscedastic/Heteroscedastic-data/Bench-", type.iter, ".RData"))
  n <- nrow(X.1)
  p.1 <- ncol(X.1)
  p.2 <- ncol(X.2)
  mu.theta <- rep(0, p.1 + p.2)
  Sigma.theta <- sigma.2.theta*diag(p.1 + p.2)
  
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
  mcmc.mu <- colMeans(mcmc.samples)
  mcmc.Sigma <- var(mcmc.samples)
  mcmc.summary <- summary(stan.res)$summary
  
  total.time <- proc.time() - start.time
  
  grid.points <- matrix(nrow = p.1 + p.2, ncol = total.grid.points)
  mcmc.values <- matrix(nrow = p.1 + p.2, ncol = total.grid.points)
  
  for (j in 1:(p.1 + p.2)) {
    grid.points[j, ] <- seq(from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                            to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                            length = total.grid.points)
    
    mcmc.values[j, ] <- demp(grid.points[j, ], obs = mcmc.samples[, j])
    
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
  mcmc.short.mu <- colMeans(mcmc.short.samples)
  mcmc.short.Sigma <- var(mcmc.short.samples)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "mcmc-short",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - demp(grid.points[j, ], 
                                                                                                              obs = mcmc.short.samples[, j])))/2,
                                                 diff_mu = (mcmc.short.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                                 diff_sigma = (sqrt(mcmc.short.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
  }
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "mcmc-short",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### EP
  
  start.time <- proc.time()
  
  ep.res <- ep.approx(X.1, X.2, y, mu.theta, Sigma.theta, 
                      eta = 0.5, alpha = 1, Q.star.init = 0.01*diag(2), r.star.init = rep(0, 2), offset = 0,
                      min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                      abs.thresh = 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "ep",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - dnorm(grid.points[j, ], 
                                                                                                               mean = ep.mu[j], 
                                                                                                               sd = sqrt(ep.Sigma[j, j]))))/2,
                                                 diff_mu = (ep.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                                 diff_sigma = (sqrt(ep.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
  }
  
  for (repetition in 1:match.reps) {
    index <- (nrow(mcmc.samples) - match.size*repetition + 1):(nrow(mcmc.samples) - match.size*(repetition - 1))
    
    bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                 method = "ep",
                                                 repetition = repetition,
                                                 match_pairs = nbp.match.pairs(rmvnorm(match.size, ep.mu, ep.Sigma),
                                                                               unname(mcmc.samples[index, ])))
  }
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "ep",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### Laplace
  
  start.time <- proc.time()
  
  laplace.res <- laplace.approx(X.1, X.2, y, mu.theta, Sigma.theta, lambda.init = 0.5, maxit = 50000)
  laplace.mu <- laplace.res$mu
  laplace.Sigma <- laplace.res$Sigma
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p.1 + p.2)) {
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "laplace",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - dnorm(grid.points[j, ], 
                                                                                                               mean = laplace.mu[j], 
                                                                                                               sd = sqrt(laplace.Sigma[j, j]))))/2,
                                                 diff_mu = (laplace.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                                 diff_sigma = (sqrt(laplace.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
  }
  
  for (repetition in 1:match.reps) {
    index <- (nrow(mcmc.samples) - match.size*repetition + 1):(nrow(mcmc.samples) - match.size*(repetition - 1))
    
    bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                 method = "laplace",
                                                 repetition = repetition,
                                                 match_pairs = nbp.match.pairs(rmvnorm(match.size, laplace.mu, laplace.Sigma),
                                                                               unname(mcmc.samples[index, ])))
  }
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "laplace",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### Means and covariances
  
  bench.res.list[[paste0("b", type.iter)]] <- list(mcmc.mu = as.vector(mcmc.mu), mcmc.Sigma = mcmc.Sigma,
                                                   mcmc.short.mu = as.vector(mcmc.short.mu), mcmc.short.Sigma = mcmc.short.Sigma,
                                                   ep.mu = as.vector(ep.mu), ep.Sigma = ep.Sigma,
                                                   laplace.mu = as.vector(laplace.mu), laplace.Sigma = laplace.Sigma)
}

save(sim.res.df.1, sim.res.df.2, sim.res.df.3, sim.res.df.4, sim.res.list,
     bench.res.df.1, bench.res.df.2, bench.res.df.3, bench.res.df.4, bench.res.list,
     file = "Heteroscedastic/Heteroscedastic-results.RData")
