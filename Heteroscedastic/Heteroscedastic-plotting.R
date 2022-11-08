
# Plots for heteroscedastic linear regression

source("EP-general-auxiliaries.R")
source("Heteroscedastic/Heteroscedastic-auxiliaries.R")

## Simulations

sim.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-Laplace.RData"))
  laplace.df <- results.df %>% mutate(method = "laplace")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, laplace.df, ep.df) %>% mutate(sim = type.iter)
  sim.plot.df <- rbind(sim.plot.df, combined.df)
}

sim.plot.df <- sim.plot.df %>% 
  mutate(sigma = sqrt(sigma_2)) %>% 
  select(-sigma_2) %>% 
  pivot_wider(names_from = method, values_from = c(mu, sigma)) %>% 
  mutate(mu_laplace = (mu_laplace - mu_mcmc)/sigma_mcmc,
         mu_ep = (mu_ep - mu_mcmc)/sigma_mcmc,
         sigma_laplace = (sigma_laplace - sigma_mcmc)/sigma_mcmc,
         sigma_ep = (sigma_ep - sigma_mcmc)/sigma_mcmc) %>% 
  select(sim, iteration, j, mu_laplace, mu_ep, sigma_laplace, sigma_ep) %>% 
  pivot_longer(cols = c(mu_laplace, mu_ep, sigma_laplace, sigma_ep)) %>% 
  separate(col = name, into = c("statistic", "method")) %>% 
  pivot_wider(names_from = statistic, values_from = value)

### Regression parameters

sim.plot.1 <- sim.plot.df %>% 
  filter(j <= as.numeric(map(sim, ~sim.settings[[.]][["p.1"]]))) %>% 
  ggplot() +
  aes(x = iteration, y = mu, colour = method) +
  geom_jitter(width = 0.3, height = 0, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~sim, nrow = num.each.type, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                      "2" = "n = 200, p.1 = 20, p.2 = 20",
                                      "3" = "n = 200, p.1 = 10, p.2 = 40"))) +
  labs(x = "Iteration",
       y = "Standardised difference in mu",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-1.png", sim.plot.1, width = 8, height = 5)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-1.png")

sim.plot.2 <- sim.plot.df %>% 
  filter(j <= as.numeric(map(sim, ~sim.settings[[.]][["p.1"]]))) %>% 
  ggplot() +
  aes(x = iteration, y = sigma, colour = method) +
  geom_jitter(width = 0.3, height = 0, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~sim, nrow = num.each.type, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                      "2" = "n = 200, p.1 = 20, p.2 = 20",
                                      "3" = "n = 200, p.1 = 10, p.2 = 40"))) +
  labs(x = "Iteration",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-2.png", sim.plot.2, width = 8, height = 5)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-2.png")

### Variance regression parameters

sim.plot.3 <- sim.plot.df %>% 
  filter(j > as.numeric(map(sim, ~sim.settings[[.]][["p.1"]]))) %>% 
  ggplot() +
  aes(x = iteration, y = mu, colour = method) +
  geom_jitter(width = 0.3, height = 0, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~sim, nrow = num.each.type, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                      "2" = "n = 200, p.1 = 20, p.2 = 20",
                                      "3" = "n = 200, p.1 = 10, p.2 = 40"))) +
  labs(x = "Iteration",
       y = "Standardised difference in mu",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-3.png", sim.plot.3, width = 8, height = 5)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-3.png")

sim.plot.4 <- sim.plot.df %>% 
  filter(j > as.numeric(map(sim, ~sim.settings[[.]][["p.1"]]))) %>% 
  ggplot() +
  aes(x = iteration, y = sigma, colour = method) +
  geom_jitter(width = 0.3, height = 0, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~sim, nrow = num.each.type, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                      "2" = "n = 200, p.1 = 20, p.2 = 20",
                                      "3" = "n = 200, p.1 = 10, p.2 = 40"))) +
  labs(x = "Iteration",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-4.png", sim.plot.4, width = 8, height = 5)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-4.png")

## Benchmarks

bench.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-Laplace.RData"))
  laplace.df <- results.df %>% mutate(method = "laplace")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, laplace.df, ep.df) %>% mutate(bench = type.iter)
  bench.plot.df <- rbind(bench.plot.df, combined.df)
}

bench.plot.df <- bench.plot.df %>% 
  mutate(sigma = sqrt(sigma_2)) %>% 
  select(-sigma_2) %>% 
  pivot_wider(names_from = method, values_from = c(mu, sigma)) %>% 
  mutate(mu_laplace = (mu_laplace - mu_mcmc)/sigma_mcmc,
         mu_ep = (mu_ep - mu_mcmc)/sigma_mcmc,
         sigma_laplace = (sigma_laplace - sigma_mcmc)/sigma_mcmc,
         sigma_ep = (sigma_ep - sigma_mcmc)/sigma_mcmc) %>% 
  select(bench, j, mu_laplace, mu_ep, sigma_laplace, sigma_ep) %>% 
  pivot_longer(cols = c(mu_laplace, mu_ep, sigma_laplace, sigma_ep)) %>% 
  separate(col = name, into = c("statistic", "method")) %>% 
  pivot_wider(names_from = statistic, values_from = value)

### Regression parameters

bench.plot.1 <- bench.plot.df %>% 
  filter(j <= as.numeric(map(bench, ~bench.settings[[.]][["p.1"]]))) %>% 
  ggplot() +
  aes(x = mu, y = sigma, colour = method) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~bench, nrow = 1, scales = "free",
             labeller = as_labeller(c("1" = "Food (n = 40, p.1 = 2, p.2 = 2)",
                                      "2" = "Salary (n = 725, p.1 = 6, p.2 = 2)",
                                      "3" = "Sniffer (n = 125, p.1 = 5, p.2 = 5)"))) +
  labs(x = "Standardised difference in mu",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-bench-plot-1.png", bench.plot.1, width = 8, height = 3)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-bench-plot-1.png")

### Variance regression parameters

bench.plot.2 <- bench.plot.df %>% 
  filter(j > as.numeric(map(bench, ~bench.settings[[.]][["p.1"]]))) %>% 
  ggplot() +
  aes(x = mu, y = sigma, colour = method) +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~bench, nrow = 1, scales = "free",
             labeller = as_labeller(c("1" = "Food (n = 40, p.1 = 2, p.2 = 2)",
                                      "2" = "Salary (n = 725, p.1 = 6, p.2 = 2)",
                                      "3" = "Sniffer (n = 125, p.1 = 5, p.2 = 5)"))) +
  labs(x = "Standardised difference in mu",
       y = "Standardised difference in sigma",
       colour = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-bench-plot-2.png", bench.plot.2, width = 8, height = 3)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-bench-plot-2.png")
