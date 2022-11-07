
# Simulations/benchmarks for quantile regression (MFVB)

source("EP-general-auxiliaries.R")
source("Quantile/Quantile-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(iteration = integer(),
                           j = double(),
                           mu = double(),
                           sigma_2 = double())
  
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    load(paste0("Quantile/Quantile-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X)
    p <- ncol(X)
    mu.theta <- rep(0, p + 1)
    Sigma.theta <- sigma.2.theta*diag(p + 1)
    
    mu.beta <- mu.theta[1:p]
    Sigma.beta <- Sigma.theta[1:p, 1:p]
    mu.kappa <- mu.theta[p + 1]
    sigma.2.kappa <- Sigma.theta[p + 1, p + 1]
    
    mfvb.res <- mfvb.approx(X, y, mu.beta, Sigma.beta, mu.kappa, sigma.2.kappa, 
                            tau, maxit = 1000, tol = 1.0E-14, verbose = F)
    
    mfvb.mu <- mfvb.res$mu
    mfvb.Sigma <- mfvb.res$Sigma
    
    for (j in 1:(p + 1)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = mfvb.mu[j],
                                           sigma_2 = mfvb.Sigma[j, j])
    }
  }
  
  save(results.df, file = paste0("Quantile/Quantile-results/Sim-", type.iter, "-res-MFVB.RData"))
}

## Benchmarks

for (type.iter in 1:num.each.type) {
  results.df <- data.frame(j = double(),
                           mu = double(),
                           sigma_2 = double())
  
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Quantile/Quantile-data/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  mu.theta <- rep(0, p + 1)
  Sigma.theta <- sigma.2.theta*diag(p + 1)
  
  mu.beta <- mu.theta[1:p]
  Sigma.beta <- Sigma.theta[1:p, 1:p]
  mu.kappa <- mu.theta[p + 1]
  sigma.2.kappa <- Sigma.theta[p + 1, p + 1]
  
  mfvb.res <- mfvb.approx(X, y, mu.beta, Sigma.beta, mu.kappa, sigma.2.kappa, 
                          tau, maxit = 1000, tol = 1.0E-14, verbose = F)

  mfvb.mu <- mfvb.res$mu
  mfvb.Sigma <- mfvb.res$Sigma
  
  for (j in 1:(p + 1)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = mfvb.mu[j],
                                         sigma_2 = mfvb.Sigma[j, j])
  }
  
  save(results.df, file = paste0("Quantile/Quantile-results/Bench-", type.iter, "-res-MFVB.RData"))
}
