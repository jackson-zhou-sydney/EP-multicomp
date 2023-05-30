
# Plots and tables for lasso linear regression

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

res.files <- list.files("Lasso/Results/")

## MCMC convergence for simulations

sim.conv.files <- res.files[grep("Simulations-conv-results", res.files)]
sim.r.hat.cdf <- data.frame()

for (file in sim.conv.files) {
  load(paste0("Lasso/Results/", file))
  sim.r.hat.cdf <- rbind(sim.r.hat.cdf, sim.r.hat.df)
}

sim.r.hat.cdf %>% 
  group_by(sim, mcmc_iter, seed) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(sim, mcmc_iter) %>% 
  summarise(mean_max_r_hat = round(mean(max_r_hat), 2))

## MCMC convergence for benchmarks
