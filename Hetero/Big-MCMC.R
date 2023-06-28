
# MCMC for heteroscedastic linear regression big data example

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
set.seed(seed)

library(cmdstanr)
library(posterior)
mcmc <- cmdstan_model("Hetero/Methods/MCMC.stan", stanc_options = list("O1"))

r.hat.df <- data.frame(mcmc_iter = integer(),
                       j = integer(),
                       r_hat = double())

load("Hetero/Data/Big/Big.RData")
n <- nrow(X.1)
p.1 <- ncol(X.1)
p.2 <- ncol(X.2)

mu.theta <- rep(0, p.1 + p.2)
Sigma.theta <- diag(c(rep(sigma.2.beta.1, p.1), rep(sigma.2.beta.2, p.2)))

start.time <- proc.time()

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
                        iter_sampling = mcmc.g.iter,
                        iter_warmup = mcmc.g.warmup,
                        max_treedepth = max.tree.depth,
                        refresh = 1)

mcmc.g.draws <- stan.res$draws()
mcmc.g.samples <- as.matrix(as_draws_df(mcmc.g.draws))[, 2:(1 + p.1 + p.2)]
tail.mcmc.g.samples <- as.matrix(as_draws_df(tail(mcmc.g.draws, round(eval.size/num.cores))))[, 2:(1 + p.1 + p.2)]
mcmc.g.mu <- colMeans(mcmc.g.samples)
mcmc.g.Sigma <- var(mcmc.g.samples)

total.time <- proc.time() - start.time
mcmc.g.time <- total.time["elapsed"]

grid.points <- matrix(nrow = p.1 + p.2, ncol = total.grid.points)
mcmc.g.values <- matrix(nrow = p.1 + p.2, ncol = total.grid.points)

for (j in 1:(p.1 + p.2)) {
  density.res <- density(mcmc.g.samples[, j], bw = "SJ-ste",
                         from = mcmc.g.mu[j] - sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                         to = mcmc.g.mu[j] + sd.multiple*sqrt(mcmc.g.Sigma[j, j]),
                         n = total.grid.points)
  
  grid.points[j, ] <- density.res$x
  mcmc.g.values[j, ] <- density.res$y
}

for (j in 1:(p.1 + p.2)) {
  param_samples <- matrix(as.numeric(mcmc.g.draws[, , 1 + j]), nrow = mcmc.g.iter)
  
  for (mcmc.iter in big.test.iter) {
    r.hat.df <- r.hat.df %>% add_row(mcmc_iter = mcmc.iter,
                                     j = j,
                                     r_hat = rstan::Rhat(param_samples[1:mcmc.iter, ]))
  }
}

r.hat.df <- r.hat.df %>% group_by(mcmc_iter) %>% summarise(max_r_hat = max(r_hat))
mcmc.s.iter <- r.hat.df %>% filter(max_r_hat < r.hat.tol) %>% pull(mcmc_iter) %>% min()

mcmc.s.draws <- head(mcmc.g.draws, mcmc.s.iter)
mcmc.s.samples <- as.matrix(as_draws_df(mcmc.s.draws))[, 2:(1 + p.1 + p.2)]
tail.mcmc.s.samples <- as.matrix(as_draws_df(tail(mcmc.s.draws, round(eval.size/num.cores))))[, 2:(1 + p.1 + p.2)]
mcmc.s.mu <- colMeans(mcmc.s.samples)
mcmc.s.Sigma <- var(mcmc.s.samples)
mcmc.s.time <- (mcmc.s.iter/mcmc.g.iter)*mcmc.g.time

save(mcmc.g.mu, mcmc.g.Sigma, tail.mcmc.g.samples, mcmc.g.time, grid.points, mcmc.g.values,
     mcmc.s.mu, mcmc.s.Sigma, tail.mcmc.s.samples, mcmc.s.time, r.hat.df, mcmc.s.iter,
     file = paste0("Hetero/Results/Big-MCMC-results-", str_pad(seed, 2, pad = "0"), ".Rdata"))
