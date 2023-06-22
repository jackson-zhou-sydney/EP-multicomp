
# Convergence evaluation for quantile regression simulations

source("General-auxiliaries.R")
source("Quantile/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
mcmc.iter <- as.numeric(args[1])
mcmc.warmup <- warmup.mult*mcmc.iter
seed <- as.numeric(args[2])
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Quantile/Methods/MCMC.stan")

sim.r.hat.df <- data.frame(seed = integer(),
                           sim = integer(),
                           mcmc_iter = integer(),
                           j = integer(),
                           r_hat = double())

for (type.iter in 1:num.sim) {
  print(paste0("Current simulation: ", type.iter))
  
  load(paste0("Quantile/Data/Simulations/Sim-", type.iter, "-iter-01.RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  mu.theta <- rep(0, p + 1)
  Sigma.theta <- sigma.2.theta*diag(p + 1)
  mu.beta <- mu.theta[1:p]
  Sigma.beta <- Sigma.theta[1:p, 1:p]
  mu.kappa <- mu.theta[p + 1]
  sigma.2.kappa <- Sigma.theta[p + 1, p + 1]
  
  stan.res <- mcmc$sample(data = list(N = n,
                                      p = p,
                                      X = X,
                                      y = y,
                                      mu_theta = mu.theta,
                                      Sigma_theta = Sigma.theta,
                                      tau = tau), 
                          seed = seed, 
                          chains = num.cores, 
                          parallel_chains = num.cores,
                          iter_sampling = mcmc.iter,
                          iter_warmup = mcmc.warmup,
                          refresh = 1)
  
  mcmc.summary <- stan.res$summary()
  
  for (j in 1:(p + 1)) {
    sim.r.hat.df <- sim.r.hat.df %>% add_row(seed = seed,
                                             sim = type.iter,
                                             mcmc_iter = mcmc.iter,
                                             j = j,
                                             r_hat = mcmc.summary %>% filter(variable == paste0("theta[", j, "]")) %>% pull(rhat))
  }
}

save(sim.r.hat.df, file = paste0("Quantile/Results/Simulations-conv-results-", mcmc.iter, "-", str_pad(seed, 2, pad = "0"), ".Rdata"))
