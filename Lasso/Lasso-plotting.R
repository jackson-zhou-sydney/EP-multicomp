
# Plots for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

## Simulations

sim.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-MFVB.RData"))
  mfvb.df <- results.df %>% mutate(method = "mfvb")
  load(paste0("Lasso/Lasso-results/Sim-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, mfvb.df, ep.df) %>% mutate(sim = type.iter)
  sim.plot.df <- rbind(sim.plot.df, combined.df)
}

sim.plot.df <- sim.plot.df %>% 
  mutate(sigma = sqrt(sigma_2)) %>% 
  select(-sigma_2) %>% 
  pivot_wider(names_from = method, values_from = c(mu, sigma)) %>% 
  mutate(mu_mfvb = (mu_mfvb - mu_mcmc)/sigma_mcmc,
         mu_ep = (mu_ep - mu_mcmc)/sigma_mcmc,
         sigma_mfvb = (sigma_mfvb - sigma_mcmc)/sigma_mcmc,
         sigma_ep = (sigma_ep - sigma_mcmc)/sigma_mcmc) %>% 
  select(sim, iteration, j, mu_mfvb, mu_ep, sigma_mfvb, sigma_ep) %>% 
  pivot_longer(cols = c(mu_mfvb, mu_ep, sigma_mfvb, sigma_ep)) %>% 
  separate(col = name, into = c("statistic", "method")) %>% 
  pivot_wider(names_from = statistic, values_from = value)

### Regression parameters

sim.plot.1 <- sim.plot.df %>% 
  filter(j <= as.numeric(map(sim, ~sim.settings[[.]][["p"]]))) %>% 
  ggplot() +
  aes(x = iteration, y = mu, colour = method) +
  geom_jitter(width = 0.3, height = 0, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "mfvb" = "MFVB"),
                      values = c("ep" = "black", "mfvb" = "darkgray")) +
  facet_wrap(~sim, nrow = num.each.type, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p = 40",
                                      "2" = "n = 40, p = 40",
                                      "3" = "n = 10, p = 40"))) +
  labs(x = "Iteration",
       y = "Standardised difference in mu",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-1.png", sim.plot.1, width = 8, height = 5)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-1.png")

sim.plot.2 <- sim.plot.df %>% 
  filter(j <= as.numeric(map(sim, ~sim.settings[[.]][["p"]]))) %>% 
  ggplot() +
  aes(x = iteration, y = sigma, colour = method) +
  geom_jitter(width = 0.3, height = 0, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "mfvb" = "MFVB"),
                      values = c("ep" = "black", "mfvb" = "darkgray")) +
  facet_wrap(~sim, nrow = num.each.type, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p = 40",
                                      "2" = "n = 40, p = 40",
                                      "3" = "n = 10, p = 40"))) +
  labs(x = "Iteration",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-2.png", sim.plot.2, width = 8, height = 5)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-2.png")

### Scale parameter

sim.plot.3 <- sim.plot.df %>% 
  filter(j == as.numeric(map(sim, ~sim.settings[[.]][["p"]])) + 1) %>% 
  ggplot() +
  aes(x = mu, y = sigma, colour = method) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "mfvb" = "MFVB"),
                      values = c("ep" = "black", "mfvb" = "darkgray")) +
  facet_wrap(~sim, nrow = 1, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p = 40",
                                      "2" = "n = 40, p = 40",
                                      "3" = "n = 10, p = 40"))) +
  labs(x = "Standardised difference in mu",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-3.png", sim.plot.3, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-3.png")

## Benchmarks

bench.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-MFVB.RData"))
  mfvb.df <- results.df %>% mutate(method = "mfvb")
  load(paste0("Lasso/Lasso-results/Bench-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, mfvb.df, ep.df) %>% mutate(bench = type.iter)
  bench.plot.df <- rbind(bench.plot.df, combined.df)
}

bench.plot.df <- bench.plot.df %>% 
  mutate(sigma = sqrt(sigma_2)) %>% 
  select(-sigma_2) %>% 
  pivot_wider(names_from = method, values_from = c(mu, sigma)) %>% 
  mutate(mu_mfvb = (mu_mfvb - mu_mcmc)/sigma_mcmc,
         mu_ep = (mu_ep - mu_mcmc)/sigma_mcmc,
         sigma_mfvb = (sigma_mfvb - sigma_mcmc)/sigma_mcmc,
         sigma_ep = (sigma_ep - sigma_mcmc)/sigma_mcmc) %>% 
  select(bench, j, mu_mfvb, mu_ep, sigma_mfvb, sigma_ep) %>% 
  pivot_longer(cols = c(mu_mfvb, mu_ep, sigma_mfvb, sigma_ep)) %>% 
  separate(col = name, into = c("statistic", "method")) %>% 
  pivot_wider(names_from = statistic, values_from = value)

### Regression parameters

bench.plot.1 <- bench.plot.df %>% 
  filter(j <= as.numeric(map(bench, ~bench.settings[[.]][["p"]]))) %>% 
  ggplot() +
  aes(x = mu, y = sigma, colour = method) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "mfvb" = "MFVB"),
                      values = c("ep" = "black", "mfvb" = "darkgray")) +
  facet_wrap(~bench, nrow = 1, scales = "free",
             labeller = as_labeller(c("1" = "Diabetes (n = 442, p = 11)",
                                      "2" = "Prostate (n = 97, p = 9)",
                                      "3" = "Eye (n = 120, p = 201)"))) +
  labs(x = "Standardised difference in mu",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-bench-plot-1.png", bench.plot.1, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-bench-plot-1.png")

### Scale parameter

bench.plot.2 <- bench.plot.df %>% 
  filter(j == as.numeric(map(bench, ~bench.settings[[.]][["p"]])) + 1) %>% 
  ggplot() +
  aes(x = mu, y = sigma, colour = method) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "mfvb" = "MFVB"),
                      values = c("ep" = "black", "mfvb" = "darkgray")) +
  facet_wrap(~bench, nrow = 1, scales = "free",
             labeller = as_labeller(c("1" = "Diabetes (n = 442, p = 11)",
                                      "2" = "Prostate (n = 97, p = 9)",
                                      "3" = "Eye (n = 120, p = 201)"))) +
  labs(x = "Standardised difference in mu",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Lasso/Lasso-plots/Lasso-bench-plot-2.png", bench.plot.2, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-bench-plot-2.png")
