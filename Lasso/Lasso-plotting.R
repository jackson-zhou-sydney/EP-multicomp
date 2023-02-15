
# Plots and tables for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")
load("Lasso/Lasso-results.RData")

## Simulations

sim.plot.1 <- sim.res.df.1 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  ggplot() +
  aes(x = method, y = -log(mmd + 0.00001)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("mcmc-a" = "MCMC-A",
                              "mcmc-b" = "MCMC-B",
                              "mcmc-c" = "MCMC-C",
                              "ep" = "EP", 
                              "mfvb" = "MFVB")) +
  facet_wrap(~sim, scales = "free_y",
             labeller = labeller(sim = as_labeller(c("1" = "n = 200, p = 40",
                                                     "2" = "n = 40, p = 40",
                                                     "3" = "n = 10, p = 40")))) +
  labs(x = "Method", y = "-log(MMD + 1e-05)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-1.png", sim.plot.1, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-1.png")

sim.plot.2 <- sim.res.df.2 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  filter(lppd != -Inf) %>% 
  group_by(sim, iteration, method) %>% 
  summarise(mean_lppd = mean(lppd)) %>% 
  ggplot() +
  aes(x = method, y = mean_lppd) +
  geom_boxplot() +
  scale_x_discrete(labels = c("mcmc-a" = "MCMC-A",
                              "mcmc-b" = "MCMC-B",
                              "mcmc-c" = "MCMC-C",
                              "ep" = "EP", 
                              "mfvb" = "MFVB")) +
  facet_wrap(~sim, scales = "free_y",
             labeller = labeller(sim = as_labeller(c("1" = "n = 200, p = 40",
                                                     "2" = "n = 40, p = 40",
                                                     "3" = "n = 10, p = 40")))) +
  labs(x = "Method", y = "Mean LPPD") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-2.png", sim.plot.2, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-2.png")

sim.table.1 <- sim.res.df.3 %>% 
  mutate(method = factor(method, levels = c("mcmc", "mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  group_by(sim, method) %>% 
  summarise(mean_time = mean(time)) %>% 
  pivot_wider(names_from = "method", values_from = "mean_time")

## Benchmarks

bench.table.1 <- bench.res.df.1 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  mutate(m_log_mmd = round(-log(mmd + 0.00001), 2)) %>% 
  dplyr::select(-mmd) %>% 
  pivot_wider(names_from = "bench", values_from = "m_log_mmd")

bench.table.2 <- bench.res.df.2 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  filter(lppd != -Inf) %>% 
  group_by(bench, method) %>% 
  summarise(mean_lppd = round(mean(lppd), 2)) %>% 
  pivot_wider(names_from = "bench", values_from = "mean_lppd")

bench.table.3 <- bench.res.df.3 %>% 
  mutate(method = factor(method, levels = c("mcmc", "mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  group_by(bench, method) %>% 
  summarise(mean_time = round(mean(time), 2)) %>%
  pivot_wider(names_from = "method", values_from = "mean_time")
