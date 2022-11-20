
# Simulations/benchmarks for heteroscedastic linear regression (EP)

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
    
    ep.res <- ep.approx(X.1, X.2, y, mu.theta, Sigma.theta, 
                        eta = 0.5, alpha = 1, Q.star.init = 0.01*diag(2), r.star.init = rep(0, 2), offset = 0,
                        min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                        abs.thresh = 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    ep.samples <- rmvnorm(match.size, ep.mu, ep.Sigma)
    
    for (j in 1:(p.1 + p.2)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = ep.mu[j],
                                           sigma_2 = ep.Sigma[j, j])
    }
    
    results.samples[[iteration]] <- ep.samples
  }
  
  save(results.df, results.samples, file = paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-EP.RData"))
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
  
  ep.res <- ep.approx(X.1, X.2, y, mu.theta, Sigma.theta, 
                      eta = 0.5, alpha = 1, Q.star.init = 0.01*diag(2), r.star.init = rep(0, 2), offset = 0,
                      min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                      abs.thresh = 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  ep.samples <- rmvnorm(match.size, ep.mu, ep.Sigma)
  
  for (j in 1:(p.1 + p.2)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = ep.mu[j],
                                         sigma_2 = ep.Sigma[j, j])
  }
  
  results.samples <- ep.samples
  
  save(results.df, results.samples, file = paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-EP.RData"))
}
