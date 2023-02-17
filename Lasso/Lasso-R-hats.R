
# R hat evaluation for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

set.seed(1)

## Simulations

sim.r.hat.df <- data.frame(sim = integer(),
                           mcmc_iter = integer(),
                           repetition = integer(),
                           j = integer(),
                           r_hat = double())

for (type.iter in 1:num.each.type) {
  print(paste0("Current progress: Simulation ", type.iter))
  
  load(paste0("Lasso/Lasso-data/Sim-", type.iter, "-iter-01.RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  for (k in 1:length(mcmc.check.iter)) {
    for (repetition in 1:r.hat.reps) {
      stan.res <- stan(file = "Lasso/Lasso-model.stan",
                       data = list(N = n,
                                   p = p,
                                   X = X,
                                   y = y,
                                   mu_kappa = mu.kappa,
                                   sigma_2_kappa = sigma.2.kappa,
                                   lambda = lambda),
                       chains = mcmc.chains,
                       iter = mcmc.check.iter[k],
                       warmup = mcmc.check.warmup[k],
                       refresh = 0,
                       init = rep(0, p + 1))
      
      mcmc.summary <- summary(stan.res)$summary
      
      for (j in 1:(p + 1)) {
        sim.r.hat.df <- sim.r.hat.df %>% add_row(sim = type.iter,
                                                 mcmc_iter = mcmc.check.iter[k],
                                                 repetition = repetition,
                                                 j = j,
                                                 r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
      }
    }
  }
}

sim.r.hat.table <- sim.r.hat.df %>% 
  group_by(sim, mcmc_iter, repetition) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(sim, mcmc_iter) %>% 
  summarise(mean_max_r_hat = mean(max_r_hat))

## Benchmarks

bench.r.hat.df <- data.frame(bench = integer(),
                             mcmc_iter = integer(),
                             repetition = integer(),
                             j = integer(),
                             r_hat = double())

for (type.iter in 1:num.each.type) {
  print(paste0("Current progress: Benchmark ", type.iter))
  
  load(paste0("Lasso/Lasso-data/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  for (k in 1:length(mcmc.check.iter)) {
    for (repetition in 1:r.hat.reps) {
      stan.res <- stan(file = "Lasso/Lasso-model.stan",
                       data = list(N = n,
                                   p = p,
                                   X = X,
                                   y = y,
                                   mu_kappa = mu.kappa,
                                   sigma_2_kappa = sigma.2.kappa,
                                   lambda = lambda),
                       chains = mcmc.chains,
                       iter = mcmc.check.iter[k],
                       warmup = mcmc.check.warmup[k],
                       refresh = 0,
                       init = rep(0, p + 1))
      
      mcmc.summary <- summary(stan.res)$summary
      
      for (j in 1:(p + 1)) {
        bench.r.hat.df <- bench.r.hat.df %>% add_row(bench = type.iter,
                                                     mcmc_iter = mcmc.check.iter[k],
                                                     repetition = repetition,
                                                     j = j,
                                                     r_hat = mcmc.summary[paste0("theta[", j, "]"), "Rhat"])
      }
    }
  }
}

bench.r.hat.table <- bench.r.hat.df %>% 
  group_by(bench, mcmc_iter, repetition) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(bench, mcmc_iter) %>% 
  summarise(mean_max_r_hat = mean(max_r_hat))

save(sim.r.hat.df, bench.r.hat.df,
     file = "Lasso/Lasso-R-hats-results.RData")
