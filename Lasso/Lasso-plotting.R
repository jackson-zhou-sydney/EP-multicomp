
# Plots for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

## Simulations

sim.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  mcmc.df <- read.csv(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MCMC.csv")) %>% mutate(method = "mcmc")
  vb.df <- read.csv(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-VB.csv")) %>% mutate(method = "vb")
  ep.df <- read.csv(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-EP.csv")) %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, vb.df, ep.df) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_vb = (mu_vb - mu_mcmc)/mu_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/mu_mcmc,
           sigma2_vb = (sigma_2_vb - sigma_2_mcmc)/sigma_2_mcmc,
           sigma2_ep = (sigma_2_ep - sigma_2_mcmc)/sigma_2_mcmc) %>% 
    select(iteration, j, mu_vb, mu_ep, sigma2_vb, sigma2_ep) %>% 
    pivot_longer(cols = c(mu_vb, mu_ep, sigma2_vb, sigma2_ep)) %>% 
    separate(col = name, into = c("statistic", "method")) %>% 
    pivot_wider(names_from = statistic, values_from = value) %>% 
    group_by(iteration, method) %>% 
    summarise(mean_mu = mean(mu),
              mean_sigma2 = mean(sigma2)) %>% 
    ungroup() %>% 
    mutate(mean_mu_sqrt = sign(mean_mu)*sqrt(abs(mean_mu)),
           mean_sigma2_sqrt = sign(mean_sigma2)*sqrt(abs(mean_sigma2)),
           sim = type.iter)
  
  sim.plot.df <- rbind(sim.plot.df, combined.df)
}

sim.plot <- ggplot(data = sim.plot.df,
                   mapping = aes(x = mean_mu_sqrt, y = mean_sigma2_sqrt, shape = method, colour = method)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(name = "Method",
                     labels = c("ep" = "EP", "vb" = "VB"),
                     values = c("ep" = 19, "vb" = 17)) +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "vb" = "VB"),
                      values = c("ep" = "black", "vb" = "darkgray")) +
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
  mcmc.df <- read.csv(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MCMC.csv")) %>% mutate(method = "mcmc")
  vb.df <- read.csv(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-VB.csv")) %>% mutate(method = "vb")
  ep.df <- read.csv(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-EP.csv")) %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, vb.df, ep.df) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_vb = (mu_vb - mu_mcmc)/mu_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/mu_mcmc,
           sigma2_vb = (sigma_2_vb - sigma_2_mcmc)/sigma_2_mcmc,
           sigma2_ep = (sigma_2_ep - sigma_2_mcmc)/sigma_2_mcmc) %>% 
    select(j, mu_vb, mu_ep, sigma2_vb, sigma2_ep) %>% 
    pivot_longer(cols = c(mu_vb, mu_ep, sigma2_vb, sigma2_ep)) %>% 
    separate(col = name, into = c("statistic", "method")) %>% 
    pivot_wider(names_from = statistic, values_from = value) %>% 
    mutate(mu_sqrt = sign(mu)*sqrt(abs(mu)),
           sigma2_sqrt = sign(sigma2)*sqrt(abs(sigma2)),
           bench = type.iter)
  
  bench.plot.df <- rbind(bench.plot.df, combined.df)
}

bench.plot <- ggplot(data = bench.plot.df,
                     mapping = aes(x = mu_sqrt, y = sigma2_sqrt, shape = method, colour = method)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(name = "Method",
                     labels = c("ep" = "EP", "vb" = "VB"),
                     values = c("ep" = 19, "vb" = 17)) +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "vb" = "VB"),
                      values = c("ep" = "black", "vb" = "darkgray")) +
  facet_wrap(~bench, scales = "free",
             labeller = as_labeller(c("1" = "Diabetes (n = 442, p = 11)",
                                      "2" = "Prostate (n = 97, p = 9)",
                                      "3" = "Eye (n = 120, p = 201)"))) +
  labs(x = "Mean difference in mu",
       y = "Mean difference in sigma squared",
       shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-bench-plots.png", bench.plot, width = 8, height = 5)
plot_crop("Lasso/Lasso-plots/Lasso-bench-plots.png")
