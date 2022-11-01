
# Simulations/benchmarks for random intercept linear regression (MCMC)

source("EP-general-auxiliaries.R")
source("Random-intercept/Random-intercept-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(iteration = integer(),
                           j = double(),
                           mu = double(),
                           sigma_2 = double())
  
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    load(paste0("Random-intercept/Random-intercept-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Z)
    mu.bk <- rep(0, p + 2)
    Sigma.bk <- sigma.2.bk*diag(p + 2)
    
    stan.res <- stan(file = "Random-intercept/Random-intercept-model.stan",
                     data = list(N = n,
                                 p = p,
                                 q = q,
                                 X = X,
                                 Z = Z,
                                 y = y,
                                 mu_bk = mu.bk,
                                 Sigma_bk = Sigma.bk),
                     chains = mcmc.chains,
                     iter = mcmc.iter,
                     warmup = mcmc.warmup,
                     refresh = 0,
                     init = rep(0, p + q + 2)) 
    
    mcmc.samples <- rstan::extract(stan.res)$theta
    mcmc.mu <- colMeans(mcmc.samples)
    mcmc.Sigma <- var(mcmc.samples)
    
    for (j in 1:(p + q + 2)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = mcmc.mu[j],
                                           sigma_2 = mcmc.Sigma[j, j])
    }
  }
  
  save(results.df, file = paste0("Random-intercept/Random-intercept-results/Sim-", type.iter, "-res-MCMC.RData"))
}

## Benchmarks

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(j = double(),
                           mu = double(),
                           sigma_2 = double())
  
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Random-intercept/Random-intercept-data/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  mu.bk <- rep(0, p + 2)
  Sigma.bk <- sigma.2.bk*diag(p + 2)
  
  stan.res <- stan(file = "Random-intercept/Random-intercept-model.stan",
                   data = list(N = n,
                               p = p,
                               q = q,
                               X = X,
                               Z = Z,
                               y = y,
                               mu_bk = mu.bk,
                               Sigma_bk = Sigma.bk),
                   chains = mcmc.chains,
                   iter = mcmc.iter,
                   warmup = mcmc.warmup,
                   refresh = 0,
                   init = rep(0, p + q + 2)) 
  
  mcmc.samples <- rstan::extract(stan.res)$theta
  mcmc.mu <- colMeans(mcmc.samples)
  mcmc.Sigma <- var(mcmc.samples)
  
  for (j in 1:(p + q + 2)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = mcmc.mu[j],
                                         sigma_2 = mcmc.Sigma[j, j])
  }
  
  save(results.df, file = paste0("Random-intercept/Random-intercept-results/Bench-", type.iter, "-res-MCMC.RData"))
}
