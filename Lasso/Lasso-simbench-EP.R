
# Simulations/benchmarks for lasso linear regression (EP)

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
    
    ep.res <- ep.approx(X, y, mu.kappa, sigma.2.kappa, 
                        lambda, eta = 0.5, alpha = 1, Q.star = 0.01*diag(2), r.star = rep(0, 2), prec = 0,
                        min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                        abs.thresh = 0.01, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = T)
    ep.mu <- ep.res$mu
    ep.Sigma <- ep.res$Sigma
    
    for (j in 1:(p + 1)) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = ep.mu[j],
                                           sigma_2 = ep.Sigma[j, j])
    }
  }
  
  write.table(results.df, file = paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-EP.csv"), row.names = F, sep = ",")
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
  
  ep.res <- ep.approx(X, y, mu.kappa, sigma.2.kappa, 
                      lambda, eta = 0.5, alpha = 1, Q.star = 0.01*diag(2), r.star = rep(0, 2), prec = 0,
                      min.passes = 6, max.passes = 200, tol.factor = Inf, stop.factor = Inf , 
                      abs.thresh = if (type.iter == 3) 10 else 0.01, rel.thresh = 0.9, delta.limit = Inf, patience = 40, verbose = T)
  ep.mu <- ep.res$mu
  ep.Sigma <- ep.res$Sigma
  
  for (j in 1:(p + 1)) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = ep.mu[j],
                                         sigma_2 = ep.Sigma[j, j])
  }
  
  write.table(results.df, file = paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-EP.csv"), row.names = F, sep = ",")
}
