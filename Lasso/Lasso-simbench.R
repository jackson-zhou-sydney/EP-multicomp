
# Simulations/benchmarks for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

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
    
    load(paste0("Lasso/Lasso-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X)
    p <- ncol(X)
    
    ### MCMC
    
    start.time <- proc.time()
    
    stan.res <- stan(file = "Lasso/Lasso-model.stan",
                     data = list(N = n,
                                 p = p,
                                 X = X,
                                 y = y,
                                 mu_kappa = mu.kappa,
                                 sigma_2_kappa = sigma.2.kappa,
                                 lambda = lambda),
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
    
    stan.res <- stan(file = "Lasso/Lasso-model.stan",
                     data = list(N = n,
                                 p = p,
                                 X = X,
                                 y = y,
                                 mu_kappa = mu.kappa,
                                 sigma_2_kappa = sigma.2.kappa,
                                 lambda = lambda),
                     chains = mcmc.chains,
                     iter = mcmc.short.iter,
                     warmup = mcmc.short.warmup,
                     refresh = 0,
                     init = rep(0, p + 1))
    
    mcmc.short.samples <- rstan::extract(stan.res)$theta
    mcmc.short.mu <- colMeans(mcmc.short.samples)
    mcmc.short.Sigma <- var(mcmc.short.samples)
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p + 1)) {
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
    
    ep.res <- ep.approx(X, y, mu.kappa, sigma.2.kappa, 
                        lambda, eta = 0.5, alpha = 1, Q.star.init = 0.01*diag(2), r.star.init = rep(0, 2), offset = 0,
                        min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                        abs.thresh = 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p + 1)) {
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
    
    ### MFVB
    
    start.time <- proc.time()
    
    mfvb.res <- mfvb.approx(X, y, mu.kappa, sigma.2.kappa,
                            lambda, maxit = 1000, tol = 1.0E-20, verbose = F)
    mfvb.mu <- mfvb.res$mu
    mfvb.Sigma <- mfvb.res$Sigma
    
    total.time <- proc.time() - start.time
    
    for (j in 1:(p + 1)) {
      sim.res.df.1 <- sim.res.df.1 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mfvb",
                                               j = j,
                                               l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - dnorm(grid.points[j, ], 
                                                                                                             mean = mfvb.mu[j], 
                                                                                                             sd = sqrt(mfvb.Sigma[j, j]))))/2,
                                               diff_mu = (mfvb.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                               diff_sigma = (sqrt(mfvb.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
    }
    
    for (repetition in 1:match.reps) {
      index <- (nrow(mcmc.samples) - match.size*repetition + 1):(nrow(mcmc.samples) - match.size*(repetition - 1))
      
      sim.res.df.2 <- sim.res.df.2 %>% add_row(sim = type.iter,
                                               iteration = iteration,
                                               method = "mfvb",
                                               repetition = repetition,
                                               match_pairs = nbp.match.pairs(rmvnorm(match.size, mfvb.mu, mfvb.Sigma),
                                                                             unname(mcmc.samples[index, ])))
    }
    
    sim.res.df.3 <- sim.res.df.3 %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mfvb",
                                             time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
    
    ### Means and covariances
    
    sim.res.list[[paste0("s", type.iter)]][[paste0("i", iteration)]] <- list(mcmc.mu = as.vector(mcmc.mu), mcmc.Sigma = mcmc.Sigma,
                                                                             mcmc.short.mu = as.vector(mcmc.short.mu), mcmc.short.Sigma = mcmc.short.Sigma,
                                                                             ep.mu = as.vector(ep.mu), ep.Sigma = ep.Sigma,
                                                                             mfvb.mu = as.vector(mfvb.mu), mfvb.Sigma = mfvb.Sigma)
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
  
  load(paste0("Lasso/Lasso-data/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  ### MCMC
  
  start.time <- proc.time()
  
  stan.res <- stan(file = "Lasso/Lasso-model.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
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
  
  stan.res <- stan(file = "Lasso/Lasso-model.stan",
                   data = list(N = n,
                               p = p,
                               X = X,
                               y = y,
                               mu_kappa = mu.kappa,
                               sigma_2_kappa = sigma.2.kappa,
                               lambda = lambda),
                   chains = mcmc.chains,
                   iter = mcmc.short.iter,
                   warmup = mcmc.short.warmup,
                   refresh = 0,
                   init = rep(0, p + 1))
  
  mcmc.short.samples <- rstan::extract(stan.res)$theta
  mcmc.short.mu <- colMeans(mcmc.short.samples)
  mcmc.short.Sigma <- var(mcmc.short.samples)
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
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
  
  ep.res <- ep.approx(X, y, mu.kappa, sigma.2.kappa, 
                      lambda, eta = 0.5, alpha = 1, Q.star.init = 0.01*diag(2), r.star.init = rep(0, 2), offset = 0,
                      min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                      abs.thresh = if (type.iter == 3) 10 else 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
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
  
  ### MFVB
  
  start.time <- proc.time()
  
  mfvb.res <- mfvb.approx(X, y, mu.kappa, sigma.2.kappa,
                          lambda, maxit = 1000, tol = 1.0E-20, verbose = F)
  mfvb.mu <- mfvb.res$mu
  mfvb.Sigma <- mfvb.res$Sigma
  
  total.time <- proc.time() - start.time
  
  for (j in 1:(p + 1)) {
    bench.res.df.1 <- bench.res.df.1 %>% add_row(bench = type.iter,
                                                 method = "mfvb",
                                                 j = j,
                                                 l1 = 1 - trapz(grid.points[j, ], abs(mcmc.values[j, ] - dnorm(grid.points[j, ], 
                                                                                                               mean = mfvb.mu[j], 
                                                                                                               sd = sqrt(mfvb.Sigma[j, j]))))/2,
                                                 diff_mu = (mfvb.mu[j] - mcmc.mu[j])/sqrt(mcmc.Sigma[j, j]),
                                                 diff_sigma = (sqrt(mfvb.Sigma[j, j]) - sqrt(mcmc.Sigma[j, j]))/sqrt(mcmc.Sigma[j, j]))
  }
  
  for (repetition in 1:match.reps) {
    index <- (nrow(mcmc.samples) - match.size*repetition + 1):(nrow(mcmc.samples) - match.size*(repetition - 1))
    
    bench.res.df.2 <- bench.res.df.2 %>% add_row(bench = type.iter,
                                                 method = "mfvb",
                                                 repetition = repetition,
                                                 match_pairs = nbp.match.pairs(rmvnorm(match.size, mfvb.mu, mfvb.Sigma),
                                                                               unname(mcmc.samples[index, ])))
  }
  
  bench.res.df.3 <- bench.res.df.3 %>% add_row(bench = type.iter,
                                               method = "mfvb",
                                               time = sum(total.time[c(1, 2, 4, 5)], na.rm = T))
  
  ### Means and covariances
  
  bench.res.list[[paste0("b", type.iter)]] <- list(mcmc.mu = as.vector(mcmc.mu), mcmc.Sigma = mcmc.Sigma,
                                                   mcmc.short.mu = as.vector(mcmc.short.mu), mcmc.short.Sigma = mcmc.short.Sigma,
                                                   ep.mu = as.vector(ep.mu), ep.Sigma = ep.Sigma,
                                                   mfvb.mu = as.vector(mfvb.mu), mfvb.Sigma = mfvb.Sigma)
}

save(sim.res.df.1, sim.res.df.2, sim.res.df.3, sim.res.df.4, sim.res.list,
     bench.res.df.1, bench.res.df.2, bench.res.df.3, bench.res.df.4, bench.res.list,
     file = "Lasso/Lasso-results.RData")
