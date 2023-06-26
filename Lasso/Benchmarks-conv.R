
# Convergence evaluation for lasso linear regression benchmarks

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
mcmc.iter <- as.numeric(args[1])
mcmc.warmup <- warmup.mult*mcmc.iter
seed <- as.numeric(args[2])
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Lasso/Methods/MCMC.stan", stanc_options = list("O1"))

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
  
  stan.res <- mcmc$sample(data = list(N = n,
                                      p = p,
                                      X = X,
                                      y = y,
                                      mu_kappa = mu.kappa,
                                      sigma_2_kappa = sigma.2.kappa,
                                      lambda = lambda), 
                          seed = seed, 
                          chains = num.cores, 
                          parallel_chains = num.cores,
                          iter_sampling = mcmc.iter,
                          iter_warmup = mcmc.warmup,
                          max_treedepth = max.tree.depth,
                          refresh = 1)
  
  mcmc.summary <- stan.res$summary()
  
  for (j in 1:(p + 1)) {
    bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 mcmc_iter = mcmc.iter,
                                                 j = j,
                                                 r_hat = mcmc.summary %>% filter(variable == paste0("theta[", j, "]")) %>% pull(rhat))
  }
}

save(bench.r.hat.df, file = paste0("Lasso/Results/Benchmarks-conv-results-", mcmc.iter, "-", str_pad(seed, 2, pad = "0"), ".Rdata"))
