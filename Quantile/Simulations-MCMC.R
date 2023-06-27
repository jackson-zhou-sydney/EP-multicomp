
# MCMC for quantile linear regression simulations

source("General-auxiliaries.R")
source("Quantile/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
set.seed(seed)

library(cmdstanr)
library(posterior)
mcmc <- cmdstan_model("Quantile/Methods/MCMC.stan", stanc_options = list("O1"))

for (type.iter in 1:num.sim) {
  print(paste0("Current simulation: ", type.iter))
  
  for (iteration in 1:num.sim.iter) {
    print(paste0("Current iteration: ", iteration))
    
    r.hat.df <- data.frame(mcmc_iter = integer(),
                           j = integer(),
                           r_hat = double())
    
    load(paste0("Quantile/Data/Simulations/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
    n <- nrow(X)
    p <- ncol(X)
    
    mu.theta <- rep(0, p + 1)
    Sigma.theta <- diag(c(rep(sigma.2.beta, p), sigma.2.kappa))
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
                            iter_sampling = mcmc.g.iter,
                            iter_warmup = mcmc.g.warmup,
                            max_treedepth = max.tree.depth,
                            refresh = 1)
    
    mcmc.g.draws <- stan.res$draws()
    mcmc.g.samples <- as.matrix(as_draws_df(mcmc.g.draws))[, 2:(1 + p + 1)]
    tail.mcmc.g.samples <- as.matrix(as_draws_df(tail(mcmc.g.draws, round(eval.size/num.cores))))[, 2:(1 + p + 1)]
    mcmc.g.mu <- colMeans(mcmc.g.samples)
    mcmc.g.Sigma <- var(mcmc.g.samples)
    
    grid.points <- matrix(nrow = p + 1, ncol = total.grid.points)
    mcmc.g.values <- matrix(nrow = p + 1, ncol = total.grid.points)
    
    for (j in 1:(p + 1)) {
      density.res <- density(mcmc.g.samples[, j], bw = "SJ-ste",
                             from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                             to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                             n = total.grid.points)
      
      grid.points[j, ] <- density.res$x
      mcmc.g.values[j, ] <- density.res$y
    }
    
    for (j in 1:(p + 1)) {
      param_samples <- matrix(as.numeric(mcmc.g.draws[, , 1 + j]), nrow = mcmc.g.iter)
      
      for (mcmc.iter in seq(from = mcmc.min.iter, to = mcmc.g.iter, by = mcmc.step)) {
        r.hat.df <- r.hat.df %>% add_row(mcmc_iter = mcmc.iter,
                                         j = j,
                                         r_hat = rstan::Rhat(param_samples[1:mcmc.iter, ]))
      }
    }
    
    r.hat.df <- r.hat.df %>% group_by(mcmc_iter) %>% summarise(max_r_hat = max(r_hat))
    mcmc.s.iter <- r.hat.df %>% filter(max_r_hat < r.hat.tol) %>% pull(mcmc_iter) %>% min()
    
    mcmc.s.draws <- head(mcmc.g.draws, mcmc.s.iter)
    mcmc.s.samples <- as.matrix(as_draws_df(mcmc.s.draws))[, 2:(1 + p + 1)]
    tail.mcmc.s.samples <- as.matrix(as_draws_df(tail(mcmc.s.draws, round(eval.size/num.cores))))[, 2:(1 + p + 1)]
    mcmc.s.mu <- colMeans(mcmc.s.samples)
    mcmc.s.Sigma <- var(mcmc.s.samples)
    
    save(mcmc.g.mu, mcmc.g.Sigma, tail.mcmc.g.samples, grid.points, mcmc.g.values,
         mcmc.s.mu, mcmc.s.Sigma, tail.mcmc.s.samples, r.hat.df, mcmc.s.iter,
         file = paste0("Quantile/Results/Simulations-MCMC-results-", type.iter, "-", str_pad(iteration, 2, pad = "0"), "-", str_pad(seed, 2, pad = "0"), ".Rdata"))
  }
}
