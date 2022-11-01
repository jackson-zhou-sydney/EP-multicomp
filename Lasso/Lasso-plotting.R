
# Plots for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

## Simulations

sim.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Lasso/Lasso-data/Sim-", type.iter, "-iter-01.RData"))
  p <- ncol(X)
  
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MFVB.RData"))
  mfvb.df <- results.df %>% mutate(method = "mfvb")
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, mfvb.df, ep.df) %>% 
    filter(j <= p) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_mfvb = (mu_mfvb - mu_mcmc)/sigma_2_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/sigma_2_mcmc,
           sigma2_mfvb = (sigma_2_mfvb - sigma_2_mcmc)/sigma_2_mcmc,
           sigma2_ep = (sigma_2_ep - sigma_2_mcmc)/sigma_2_mcmc) %>% 
    select(iteration, j, mu_mfvb, mu_ep, sigma2_mfvb, sigma2_ep) %>% 
    pivot_longer(cols = c(mu_mfvb, mu_ep, sigma2_mfvb, sigma2_ep)) %>% 
    separate(col = name, into = c("statistic", "method")) %>% 
    pivot_wider(names_from = statistic, values_from = value) %>% 
    group_by(iteration, method) %>% 
    summarise(mean_mu = mean(mu),
              mean_sigma2 = mean(sigma2)) %>% 
    ungroup() %>% 
    mutate(sim = type.iter)
  
  sim.plot.df <- rbind(sim.plot.df, combined.df)
}

sim.plot <- ggplot(data = sim.plot.df,
                   mapping = aes(x = mean_mu, y = mean_sigma2, shape = method, colour = method)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(name = "Method",
                     labels = c("ep" = "EP", "mfvb" = "MFVB"),
                     values = c("ep" = 19, "mfvb" = 17)) +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "mfvb" = "MFVB"),
                      values = c("ep" = "black", "mfvb" = "darkgray")) +
  facet_wrap(~sim, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p = 40",
                                      "2" = "n = 40, p = 40",
                                      "3" = "n = 10, p = 40"))) +
  labs(x = "Mean difference in mu",
       y = "Mean difference in sigma squared",
       shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-sim-plots.png", sim.plot, width = 8, height = 5)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plots.png")

## Benchmarks

bench.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Lasso/Lasso-data/Bench-", type.iter, ".RData"))
  p <- ncol(X)
  
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MFVB.RData"))
  mfvb.df <- results.df %>% mutate(method = "mfvb")
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, mfvb.df, ep.df) %>% 
    filter(j <= p) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_mfvb = (mu_mfvb - mu_mcmc)/sigma_2_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/sigma_2_mcmc,
           sigma2_mfvb = (sigma_2_mfvb - sigma_2_mcmc)/sigma_2_mcmc,
           sigma2_ep = (sigma_2_ep - sigma_2_mcmc)/sigma_2_mcmc) %>% 
    select(j, mu_mfvb, mu_ep, sigma2_mfvb, sigma2_ep) %>% 
    pivot_longer(cols = c(mu_mfvb, mu_ep, sigma2_mfvb, sigma2_ep)) %>% 
    separate(col = name, into = c("statistic", "method")) %>% 
    pivot_wider(names_from = statistic, values_from = value) %>% 
    mutate(bench = type.iter)
  
  bench.plot.df <- rbind(bench.plot.df, combined.df)
}

bench.plot <- ggplot(data = bench.plot.df,
                     mapping = aes(x = mu, y = sigma2, shape = method, colour = method)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(name = "Method",
                     labels = c("ep" = "EP", "mfvb" = "MFVB"),
                     values = c("ep" = 19, "mfvb" = 17)) +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "mfvb" = "MFVB"),
                      values = c("ep" = "black", "mfvb" = "darkgray")) +
  facet_wrap(~bench, scales = "free",
             labeller = as_labeller(c("1" = "Diabetes (n = 442, p = 11)",
                                      "2" = "Prostate (n = 97, p = 9)",
                                      "3" = "Eye (n = 120, p = 201)"))) +
  labs(x = "Difference in mu",
       y = "Difference in sigma squared",
       shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-bench-plots.png", bench.plot, width = 8, height = 5)
plot_crop("Lasso/Lasso-plots/Lasso-bench-plots.png")
