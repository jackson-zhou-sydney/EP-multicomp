
# R hat evaluation for lasso linear regression benchmarks

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

set.seed(1)

bench.r.hat.df <- data.frame(bench = integer(),
                             mcmc_iter = integer(),
                             repetition = integer(),
                             j = integer(),
                             r_hat = double())

for (type.iter in 1:num.bench) {
  print(paste0("Current benchmark: ", type.iter))
  
  load(paste0("Lasso/Data/Benchmarks/Bench-", type.iter, ".RData"))
  n <- nrow(X)
  p <- ncol(X)
  
  for (k in 1:length(mcmc.check.iter)) {
    print(paste0("Current check.iter: ", mcmc.check.iter[k]))
    for (repetition in 1:r.hat.reps) {
      print(paste0("Current repetition: ", repetition))
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
  summarise(mean_max_r_hat = round(mean(max_r_hat), 2))

save(bench.r.hat.df, file = "Lasso/Results/Benchmarks-conv-results.RData")
