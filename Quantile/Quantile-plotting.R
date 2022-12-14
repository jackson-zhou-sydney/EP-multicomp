
# Plots and tables for quantile regression

source("EP-general-auxiliaries.R")
source("Quantile/Quantile-auxiliaries.R")
load("Quantile/Quantile-results.RData")

## Simulations

sim.plot.1 <- sim.res.df.1 %>% 
  mutate(method = factor(method, levels = c("mcmc-short", "ep", "mfvb"))) %>% 
  mutate(block = case_when(j <= as.numeric(map(sim, ~sim.settings[[.]][["p"]])) ~ "beta",
                           TRUE ~ "kappa")) %>% 
  group_by(sim, iteration, block, method) %>% 
  summarise(mean_l1 = mean(l1)) %>% 
  ggplot() +
  aes(x = method, y = mean_l1) +
  geom_boxplot() +
  scale_x_discrete(labels = c("mcmc-short" = "MCMC-S",
                              "ep" = "EP",
                              "mfvb" = "MFVB")) +
  facet_wrap(block ~ sim, scales = "free_y",
             labeller = labeller(block = as_labeller(c("beta" = "β", 
                                                       "kappa" = "K")),
                                 sim = as_labeller(c("1" = "Normal (n = 200, p = 40)",
                                                     "2" = "Poisson (n = 200, p = 40)",
                                                     "3" = "Binomial (n = 200, p = 40)")))) +
  labs(x = "Method", y = "Mean L1 accuracy") +
  theme_bw()

ggsave("Quantile/Quantile-plots/Quantile-sim-plot-1.png", sim.plot.1, width = 8, height = 5)
plot_crop("Quantile/Quantile-plots/Quantile-sim-plot-1.png")

sim.plot.2 <- sim.res.df.2 %>% 
  mutate(method = factor(method, levels = c("mcmc-short", "ep", "mfvb"))) %>% 
  group_by(sim, iteration, method) %>% 
  summarise(mean_match_pairs = mean(match_pairs)) %>% 
  ggplot() +
  aes(x = method, y = mean_match_pairs) +
  geom_boxplot() +
  scale_x_discrete(labels = c("mcmc-short" = "MCMC-S",
                              "ep" = "EP",
                              "mfvb" = "MFVB")) +
  facet_wrap(~sim, scales = "free_y",
             labeller = labeller(sim = as_labeller(c("1" = "Normal (n = 200, p = 40)",
                                                     "2" = "Poisson (n = 200, p = 40)",
                                                     "3" = "Binomial (n = 200, p = 40)")))) +
  labs(x = "Method", y = "Mean match pairs") +
  theme_bw()

ggsave("Quantile/Quantile-plots/Quantile-sim-plot-2.png", sim.plot.2, width = 8, height = 3)
plot_crop("Quantile/Quantile-plots/Quantile-sim-plot-2.png")

sim.table.1 <- sim.res.df.3 %>% 
  group_by(sim, method) %>% 
  summarise(mean_time = mean(time)) %>% 
  pivot_wider(names_from = "method", values_from = "mean_time")

## Benchmarks

bench.table.1 <- bench.res.df.1 %>% 
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p"]])) ~ "beta",
                           TRUE ~ "kappa")) %>% 
  group_by(bench, block, method) %>% 
  summarise(mean_l1 = mean(l1)) %>% 
  pivot_wider(names_from = "bench", values_from = "mean_l1")

bench.table.2 <- bench.res.df.2 %>% 
  group_by(bench, method) %>% 
  summarise(mean_match_pairs = mean(match_pairs)) %>% 
  pivot_wider(names_from = "bench", values_from = "mean_match_pairs")

bench.table.3 <- bench.res.df.3 %>% 
  group_by(bench, method) %>% 
  summarise(mean_time = mean(time)) %>% 
  pivot_wider(names_from = "method", values_from = "mean_time")
