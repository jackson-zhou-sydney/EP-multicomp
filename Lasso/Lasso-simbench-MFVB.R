
# Simulations/benchmarks for lasso linear regression (MFVB)

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
    
    mfvb.res <- mfvb.approx(X, y, mu.kappa, sigma.2.kappa,
                            lambda, max.iter = 1000, tol = 1.0E-10, verbose = F)

    mfvb.mu <- mfvb.res$mu
    mfvb.Sigma <- mfvb.res$Sigma
    
    for (j in 1:p) {
      results.df <- results.df %>% add_row(iteration = iteration,
                                           j = j,
                                           mu = mfvb.mu[j],
                                           sigma_2 = mfvb.Sigma[j, j])
    }
  }
  
  write.table(results.df, file = paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MFVB.csv"), row.names = F, sep = ",")
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
  
  mfvb.res <- mfvb.approx(X, y, mu.kappa, sigma.2.kappa,
                          lambda, max.iter = 1000, tol = 1.0E-10, verbose = F)

  mfvb.mu <- mfvb.res$mu
  mfvb.Sigma <- mfvb.res$Sigma
  
  for (j in 1:p) {
    results.df <- results.df %>% add_row(j = j,
                                         mu = mfvb.mu[j],
                                         sigma_2 = mfvb.Sigma[j, j])
  }
  
  write.table(results.df, file = paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MFVB.csv"), row.names = F, sep = ",")
}
