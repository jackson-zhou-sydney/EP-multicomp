
# Simulations/benchmarks for quantile regression (MCMC)

source("EP-general-auxiliaries.R")
source("Quantile/Quantile-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(iteration = integer(),
                           j = double(),
                           mu = double(),
                           sigma_2 = double(),
                           r_hat = double())
  
  times.df <- data.frame(iteration = integer(),
                         time = double())
  
  results.values <- vector(mode = "list", length = num.sim)
  results.samples <- vector(mode = "list", length = num.sim)
  
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    load(paste0("Quantile/Quantile-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X)
    p <- ncol(X)
    mu.theta <- rep(0, p + 1)
    Sigma.theta <- sigma.2.theta*diag(p + 1)
    
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
    
    mcmc.values <- matrix(nrow = p + 1, ncol = total.grid.points)
    
    for (j in 1:(p + 1)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = mcmc.mu[j],
                                           sigma_2 = mcmc.Sigma[j, j],
                                           r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
      
      mcmc.values[j, ] <- demp(seq(from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                                   to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                                   length = total.grid.points), 
                               obs = mcmc.samples[, j])
    }
    
    times.df <- times.df %>% add_row(iteration = iteration,
                                     time = sum(total.time[c(1, 2, 4, 5)], na.rm = TRUE))
    
    results.values[[iteration]] <- mcmc.values
    results.samples[[iteration]] <- unname(tail(mcmc.samples, match.size))
  }
  
  save(results.df, times.df, results.samples, results.values,
       file = paste0("Quantile/Quantile-results/Sim-", type.iter, "-res-MCMC.RData"))
}

## Benchmarks

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(j = double(),
                           mu = double(),
                           sigma_2 = double(),
                           r_hat = double())
  
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Quantile/Quantile-data/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  mu.theta <- rep(0, p + 1)
  Sigma.theta <- sigma.2.theta*diag(p + 1)
  
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
  
  mcmc.values <- matrix(nrow = p + 1, ncol = total.grid.points)
  
  for (j in 1:(p + 1)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = mcmc.mu[j],
                                         sigma_2 = mcmc.Sigma[j, j],
                                         r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
    
    mcmc.values[j, ] <- demp(seq(from = mcmc.mu[j] - sd.multiple*sqrt(mcmc.Sigma[j, j]),
                                 to = mcmc.mu[j] + sd.multiple*sqrt(mcmc.Sigma[j, j]),
                                 length = total.grid.points), 
                             obs = mcmc.samples[, j])
  }
  
  times.df <- data.frame(time = sum(total.time[c(1, 2, 4, 5)], na.rm = TRUE))
  results.values <- mcmc.values
  results.samples <- unname(tail(mcmc.samples, match.size))
  
  save(results.df, times.df, results.samples, results.values,
       file = paste0("Quantile/Quantile-results/Bench-", type.iter, "-res-MCMC.RData"))
}
