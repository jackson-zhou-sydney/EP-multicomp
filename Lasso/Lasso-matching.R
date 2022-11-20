
# Matching for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

## Simulations

sim.match.df <- data.frame(sim = double(),
                           iteration = double(),
                           method = character(),
                           match_pairs = double())

for (type.iter in 1:num.each.type) {
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MCMC.RData"))
  print(paste0("Maximum R hat for simulation ", type.iter, " is ", max(results.df$r_hat)))
  mcmc.samples <- results.samples
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MFVB.RData"))
  mfvb.samples <- results.samples
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-EP.RData"))
  ep.samples <- results.samples
  
  for (iteration in 1:num.sim) {
    print(paste0("Current progress: Simulation ", type.iter, ", iteration ", iteration, " of ", num.sim))
    
    sim.match.df <- sim.match.df %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "mfvb",
                                             match_pairs = nbp.match.pairs(mfvb.samples[[iteration]],
                                                                           mcmc.samples[[iteration]]))
    
    sim.match.df <- sim.match.df %>% add_row(sim = type.iter,
                                             iteration = iteration,
                                             method = "ep",
                                             match_pairs = nbp.match.pairs(ep.samples[[iteration]],
                                                                           mcmc.samples[[iteration]]))
  }
}

sim.mplot <-  sim.match.df %>% 
  ggplot() +
  aes(x = method, y = match_pairs, fill = method) +
  geom_boxplot() +
  scale_fill_manual(name = "Method",
                    labels = c("ep" = "EP", "mfvb" = "MFVB"),
                    values = c("ep" = "gray", "mfvb" = "white")) +
  scale_x_discrete(labels = c("ep" = "EP", "mfvb" = "MFVB")) +
  facet_wrap(~sim, nrow = 1,
             labeller = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                      "2" = "n = 200, p.1 = 20, p.2 = 20",
                                      "3" = "n = 200, p.1 = 10, p.2 = 40"))) +
  labs(x = "Method", y = "NBP matching pairs", fill = "Method") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("Lasso/Lasso-plots/Lasso-sim-mplot.png", sim.mplot, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-sim-mplot.png")

## Benchmarks

bench.match.df <- data.frame(bench = double(),
                             method = character(),
                             match_pairs = double())

for (type.iter in 1:num.each.type) {
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MCMC.RData"))
  print(paste0("Maximum R hat for benchmark ", type.iter, " is ", max(results.df$r_hat)))
  mcmc.samples <- results.samples
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MFVB.RData"))
  mfvb.samples <- results.samples
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-EP.RData"))
  ep.samples <- results.samples
  
  print(paste0("Current progress: Benchmark ", type.iter))
  
  bench.match.df <- bench.match.df %>% add_row(bench = type.iter,
                                               method = "mfvb",
                                               match_pairs = nbp.match.pairs(mfvb.samples,
                                                                             mcmc.samples))
  
  bench.match.df <- bench.match.df %>% add_row(bench = type.iter,
                                               method = "ep",
                                               match_pairs = nbp.match.pairs(ep.samples,
                                                                             mcmc.samples))
}
