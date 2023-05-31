
# Convergence evaluation for lasso linear regression simulations

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Lasso/Models/MCMC.stan")

sim.r.hat.df <- data.frame(seed = integer(),
                           sim = integer(),
                           mcmc_iter = integer(),
                           j = integer(),
                           r_hat = double())

for (type.iter in 1:num.sim) {
  print(paste0("Current simulation: ", type.iter))
  
  load(paste0("Lasso/Data/Simulations/Sim-", type.iter, "-iter-01.RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  for (k in 1:length(mcmc.check.iter)) {
    print(paste0("Current check.iter: ", mcmc.check.iter[k]))
    
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
                            iter_sampling = mcmc.check.iter[k],
                            iter_warmup = mcmc.check.warmup[k],
                            refresh = 100)
    
    mcmc.summary <- stan.res$summary()
    
    for (j in 1:(p + 1)) {
      sim.r.hat.df <- sim.r.hat.df %>% add_row(seed = seed,
                                               sim = type.iter,
                                               mcmc_iter = mcmc.check.iter[k],
                                               j = j,
                                               r_hat = mcmc.summary %>% filter(variable == paste0("theta[", j, "]")) %>% pull(rhat))
    }
  }
}

save(sim.r.hat.df, file = paste0("Lasso/Results/Simulations-conv-results-", str_pad(seed, 2, pad = "0"), ".Rdata"))
