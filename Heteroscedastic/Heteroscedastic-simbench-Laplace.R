
# Simulations/benchmarks for heteroscedastic linear regression (Laplace)

source("EP-general-auxiliaries.R")
source("Heteroscedastic/Heteroscedastic-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(iteration = integer(),
                           j = double(),
                           mu = double(),
                           sigma_2 = double())
  results.samples <- vector(mode = "list", length = num.sim)
  
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    load(paste0("Heteroscedastic/Heteroscedastic-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X.1)
    p.1 <- ncol(X.1)
    p.2 <- ncol(X.2)
    mu.theta <- rep(0, p.1 + p.2)
    Sigma.theta <- sigma.2.theta*diag(p.1 + p.2)
    
    laplace.res <- laplace.approx(X.1, X.2, y, mu.theta, Sigma.theta, lambda.init = 0.5, maxit = 50000)
    laplace.mu <- laplace.res$mu
    laplace.Sigma <- laplace.res$Sigma
    laplace.samples <- rmvnorm(match.size, laplace.mu, laplace.Sigma)
    
    for (j in 1:(p.1 + p.2)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = laplace.mu[j],
                                           sigma_2 = laplace.Sigma[j, j])
    }
    
    results.samples[[iteration]] <- laplace.samples
  }
  
  save(results.df, results.samples, file = paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-Laplace.RData"))
}

## Benchmarks

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(j = double(),
                           mu = double(),
                           sigma_2 = double())
  
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Heteroscedastic/Heteroscedastic-data/Bench-", type.iter, ".RData"))
  n <- nrow(X.1)
  p.1 <- ncol(X.1)
  p.2 <- ncol(X.2)
  mu.theta <- rep(0, p.1 + p.2)
  Sigma.theta <- sigma.2.theta*diag(p.1 + p.2)
  
  laplace.res <- laplace.approx(X.1, X.2, y, mu.theta, Sigma.theta, lambda.init = 0.5, maxit = 50000)
  laplace.mu <- laplace.res$mu
  laplace.Sigma <- laplace.res$Sigma
  laplace.samples <- rmvnorm(match.size, laplace.mu, laplace.Sigma)
  
  for (j in 1:(p.1 + p.2)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = laplace.mu[j],
                                         sigma_2 = laplace.Sigma[j, j])
  }
  
  results.samples <- laplace.samples
  
  save(results.df, results.samples, file = paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-Laplace.RData"))
}
