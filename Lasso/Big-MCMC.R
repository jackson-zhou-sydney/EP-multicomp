
# MCMC for lasso linear regression big data example

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
set.seed(seed)

library(cmdstanr)
library(posterior)
mcmc <- cmdstan_model("Lasso/Methods/MCMC.stan", stanc_options = list("O1"))

r.hat.df <- data.frame(mcmc_iter = integer(),
                       j = integer(),
                       r_hat = double())

load("Lasso/Data/Big/Big.RData")
n <- nrow(X)
p <- ncol(X)

start.time <- proc.time()

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
                        iter_sampling = big.mcmc.g.iter,
                        iter_warmup = big.mcmc.g.warmup,
                        max_treedepth = max.tree.depth,
                        init = 0,
                        refresh = 1)

mcmc.g.draws <- stan.res$draws()
mcmc.g.samples <- as.matrix(as_draws_df(mcmc.g.draws))[, 2:(1 + p + 1)]
tail.mcmc.g.samples <- as.matrix(as_draws_df(tail(mcmc.g.draws, round(eval.size/num.cores))))[, 2:(1 + p + 1)]
mcmc.g.mu <- colMeans(mcmc.g.samples)
mcmc.g.Sigma <- var(mcmc.g.samples)

total.time <- proc.time() - start.time
mcmc.g.time <- total.time["elapsed"]

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
  param_samples <- matrix(as.numeric(mcmc.g.draws[, , 1 + j]), nrow = big.mcmc.g.iter)
  
  for (mcmc.iter in big.test.iter) {
    r.hat.df <- r.hat.df %>% add_row(mcmc_iter = mcmc.iter,
                                     j = j,
                                     r_hat = rstan::Rhat(param_samples[1:mcmc.iter, ]))
  }
}

r.hat.df <- r.hat.df %>% group_by(mcmc_iter) %>% summarise(max_r_hat = max(r_hat))
iters <- r.hat.df %>% filter(max_r_hat > r.hat.tol) %>% pull(mcmc_iter)
ind <- if(length(iters) != 0) which(big.test.iter == max(iters)) + 1 else 1
mcmc.s.iter <- big.test.iter[min(ind, length(big.test.iter))]

mcmc.s.draws <- head(mcmc.g.draws, mcmc.s.iter)
mcmc.s.samples <- as.matrix(as_draws_df(mcmc.s.draws))[, 2:(1 + p + 1)]
tail.mcmc.s.samples <- as.matrix(as_draws_df(tail(mcmc.s.draws, round(eval.size/num.cores))))[, 2:(1 + p + 1)]
mcmc.s.mu <- colMeans(mcmc.s.samples)
mcmc.s.Sigma <- var(mcmc.s.samples)
mcmc.s.time <- (mcmc.s.iter/big.mcmc.g.iter)*mcmc.g.time

mcmc.s.grid.points <- matrix(nrow = p + 1, ncol = total.grid.points)
mcmc.s.values <- matrix(nrow = p + 1, ncol = total.grid.points)

for (j in 1:(p + 1)) {
  density.res <- density(mcmc.s.samples[, j], bw = "SJ-ste",
                         from = mcmc.s.mu[j] - sd.multiple*sqrt(mcmc.s.Sigma[j, j]),
                         to = mcmc.s.mu[j] + sd.multiple*sqrt(mcmc.s.Sigma[j, j]),
                         n = total.grid.points)
  
  mcmc.s.grid.points[j, ] <- density.res$x
  mcmc.s.values[j, ] <- density.res$y
}

save(mcmc.g.mu, mcmc.g.Sigma, tail.mcmc.g.samples, mcmc.g.time, grid.points, mcmc.g.values,
     mcmc.s.mu, mcmc.s.Sigma, tail.mcmc.s.samples, mcmc.s.time, mcmc.s.grid.points, mcmc.s.values, r.hat.df, mcmc.s.iter,
     file = paste0("Lasso/Results/Big-MCMC-results-", str_pad(seed, 2, pad = "0"), ".Rdata"))
