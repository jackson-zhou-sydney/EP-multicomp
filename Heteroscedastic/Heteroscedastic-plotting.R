
# Plots and tables for heteroscedastic linear regression

source("EP-general-auxiliaries.R")
source("Heteroscedastic/Heteroscedastic-auxiliaries.R")
load("Heteroscedastic/Heteroscedastic-results.RData")

## Simulations

sim.plot.1 <- sim.res.df.1 %>% 
  mutate(block = case_when(j <= as.numeric(map(sim, ~sim.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2")) %>% 
  group_by(sim, iteration, block, method) %>% 
  summarise(mean_l1 = mean(l1)) %>% 
  ggplot() +
  aes(x = method, y = mean_l1) +
  geom_boxplot() +
  scale_x_discrete(labels = c("ep" = "EP", "laplace" = "Laplace")) +
  facet_grid(rows = vars(block), cols = vars(sim),
             labeller = labeller(block = as_labeller(c("beta_1" = "β.1", 
                                                       "beta_2" = "β.2")),
                                 sim = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                                     "2" = "n = 200, p.1 = 20, p.2 = 20",
                                                     "3" = "n = 200, p.1 = 10, p.2 = 40")))) +
  labs(x = "Method", y = "Mean L1 accuracy") +
  theme_bw()

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-1.png", sim.plot.1, width = 8, height = 5)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-1.png")

sim.plot.2 <- sim.res.df.2 %>% 
  group_by(sim, iteration, method) %>% 
  summarise(mean_match_pairs = mean(match_pairs)) %>% 
  ggplot() +
  aes(x = method, y = mean_match_pairs) +
  geom_boxplot() +
  scale_x_discrete(labels = c("ep" = "EP", "laplace" = "Laplace")) +
  facet_grid(cols = vars(sim),
             labeller = labeller(sim = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                                     "2" = "n = 200, p.1 = 20, p.2 = 20",
                                                     "3" = "n = 200, p.1 = 10, p.2 = 40")))) +
  labs(x = "Method", y = "Mean match pairs") +
  theme_bw()

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-2.png", sim.plot.2, width = 8, height = 3)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plot-2.png")

sim.table.1 <- sim.res.df.3 %>% 
  group_by(sim, method) %>% 
  summarise(mean_time = mean(time)) %>% 
  pivot_wider(names_from = "method", values_from = "mean_time")

## Benchmarks

bench.table.1 <- bench.res.df.1 %>% 
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2")) %>% 
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
