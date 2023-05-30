
# R hat evaluation for lasso linear regression benchmarks

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
set.seed(seed)

bench.r.hat.df <- data.frame(seed = integer(),
                             bench = integer(),
                             mcmc_iter = integer(),
                             j = integer(),
                             r_hat = double())

for (type.iter in 1:num.bench) {
  print(paste0("Current benchmark: ", type.iter))
  
  load(paste0("Lasso/Data/Benchmarks/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  for (k in 1:length(mcmc.check.iter)) {
    print(paste0("Current check.iter: ", mcmc.check.iter[k]))
    
    stan.res <- stan(file = "Lasso/Models/MCMC.stan",
                     data = list(N = n,
                                 p = p,
                                 X = X,
                                 y = y,
                                 mu_kappa = mu.kappa,
                                 sigma_2_kappa = sigma.2.kappa,
                                 lambda = lambda),
                     chains = num.cores,
                     cores = num.cores,
                     iter = mcmc.check.iter[k],
                     warmup = mcmc.check.warmup[k],
                     refresh = 0,
                     init = rep(0, p + 1))
    
    mcmc.summary <- summary(stan.res)$summary
    
    for (j in 1:(p + 1)) {
      bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                   bench = type.iter,
                                                   mcmc_iter = mcmc.check.iter[k],
                                                   j = j,
                                                   r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
    }
  }
}

save(bench.r.hat.df, file = paste0("Lasso/Results/Benchmarks-conv-results-", str_pad(seed, 2, pad = "0"), ".Rdata"))
