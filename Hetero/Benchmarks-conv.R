
# Convergence evaluation for heteroscedastic linear regression benchmarks

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
mcmc.iter <- as.numeric(args[1])
mcmc.warmup <- warmup.mult*mcmc.iter
seed <- as.numeric(args[2])
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Hetero/Methods/MCMC.stan", stanc_options = list("O1"))

bench.r.hat.df <- data.frame(seed = integer(),
                             bench = integer(),
                             mcmc_iter = integer(),
                             j = integer(),
                             r_hat = double())

for (type.iter in 1:num.bench) {
  print(paste0("Current benchmark: ", type.iter))
  
  load(paste0("Hetero/Data/Benchmarks/Bench-", type.iter, ".RData"))
  n <- nrow(X.1)
  p.1 <- ncol(X.1)
  p.2 <- ncol(X.2)
  
  mu.theta <- rep(0, p.1 + p.2)
  Sigma.theta <- diag(c(rep(sigma.2.beta.1, p.1), rep(sigma.2.beta.2, p.2)))
  
  stan.res <- mcmc$sample(data = list(N = n,
                                      p_1 = p.1,
                                      p_2 = p.2,
                                      X_1 = X.1,
                                      X_2 = X.2,
                                      y = y,
                                      mu_theta = mu.theta,
                                      Sigma_theta = Sigma.theta),
                          seed = seed, 
                          chains = num.cores, 
                          parallel_chains = num.cores,
                          iter_sampling = mcmc.iter,
                          iter_warmup = mcmc.warmup,
                          max_treedepth = max.tree.depth,
                          init = 0,
                          refresh = 1)
  
  mcmc.summary <- stan.res$summary()
  
  for (j in 1:(p.1 + p.2)) {
    bench.r.hat.df <- bench.r.hat.df %>% add_row(seed = seed,
                                                 bench = type.iter,
                                                 mcmc_iter = mcmc.iter,
                                                 j = j,
                                                 r_hat = mcmc.summary %>% filter(variable == paste0("theta[", j, "]")) %>% pull(rhat))
  }
}

save(bench.r.hat.df, file = paste0("Hetero/Results/Benchmarks-conv-results-", mcmc.iter, "-", str_pad(seed, 2, pad = "0"), ".Rdata"))
