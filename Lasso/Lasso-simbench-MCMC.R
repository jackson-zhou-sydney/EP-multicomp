
# Simulations/benchmarks for lasso linear regression (MCMC)

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(iteration = integer(),
                           j = double(),
                           mu = double(),
                           sigma_2 = double())
  
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    df <- read.csv(paste0("Lasso/Lasso-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".csv"), header = F)
    n <- nrow(df)
    p <- ncol(df) - 1
    X <- unname(as.matrix(df[, 1:p]))
    y <- as.numeric(df[, p + 1])
    
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
                     init = "random")
    
    mcmc.samples <- rstan::extract(stan.res)$theta
    mcmc.mu <- colMeans(mcmc.samples)
    mcmc.Sigma <- var(mcmc.samples)
    
    for (j in 1:(p + 1)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = mcmc.mu[j],
                                           sigma_2 = mcmc.Sigma[j, j])
    }
  }
  
  write.table(results.df, file = paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MCMC.csv"), row.names = F, sep = ",")
}

## Benchmarks

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(j = double(),
                           mu = double(),
                           sigma_2 = double())
  
  print(paste0("Current progress: Benchmark ", type.iter))
  
  df <- read.csv(paste0("Lasso/Lasso-data/Bench-", type.iter, ".csv"), header = F)
  n <- nrow(df)
  p <- ncol(df) - 1
  X <- unname(as.matrix(df[, 1:p]))
  y <- as.numeric(df[, p + 1])
  
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
                   init = "random")
  
  mcmc.samples <- rstan::extract(stan.res)$theta
  mcmc.mu <- colMeans(mcmc.samples)
  mcmc.Sigma <- var(mcmc.samples)
  
  for (j in 1:(p + 1)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = mcmc.mu[j],
                                         sigma_2 = mcmc.Sigma[j, j])
  }
  
  write.table(results.df, file = paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MCMC.csv"), row.names = F, sep = ",")
}
