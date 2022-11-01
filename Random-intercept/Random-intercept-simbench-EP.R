
# Simulations/benchmarks for random intercept linear regression (EP)

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
    
    ep.res <- ep.approx(X, Z, y, mu.bk, Sigma.bk,
                        eta = 0.5, alpha = 1, Q.star = 0.01*diag(2), r.star = rep(0, 2), prec = 0,
                        min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf, 
                        abs.thresh = 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    
    for (j in 1:(p + q + 2)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = ep.mu[j],
                                           sigma_2 = ep.Sigma[j, j])
    }
  }
  
  save(results.df, file = paste0("Random-intercept/Random-intercept-results/Sim-", type.iter, "-res-EP.RData"))
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
  
  ep.res <- ep.approx(X, Z, y, mu.bk, Sigma.bk,
                      eta = 0.5, alpha = 1, Q.star = 0.01*diag(2), r.star = rep(0, 2), prec = 0,
                      min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf, 
                      abs.thresh = 0.1, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = F)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  
  for (j in 1:(p + q + 2)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = ep.mu[j],
                                         sigma_2 = ep.Sigma[j, j])
  }
  
  save(results.df, file = paste0("Random-intercept/Random-intercept-results/Bench-", type.iter, "-res-EP.RData"))
}
