
# Plots and tables for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")
load("Lasso/Lasso-results.RData")

## Simulations

sim.plot.1 <- sim.res.df.1 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  mutate(block = case_when(j <= as.numeric(map(sim, ~sim.settings[[.]][["p"]])) ~ "beta",
                           TRUE ~ "kappa")) %>% 
  group_by(sim, iteration, block, method) %>% 
  summarise(mean_l1 = mean(l1)) %>% 
  ggplot() +
  aes(x = method, y = mean_l1) +
  geom_boxplot() +
  scale_x_discrete(labels = c("mcmc-a" = "MCMC-A",
                              "mcmc-b" = "MCMC-B",
                              "mcmc-c" = "MCMC-C",
                              "ep" = "EP", 
                              "mfvb" = "MFVB")) +
  facet_wrap(block ~ sim, scales = "free_y",
             labeller = labeller(block = as_labeller(c("beta" = "β", 
                                                       "kappa" = "K")),
                                 sim = as_labeller(c("1" = "n = 200, p = 40",
                                                     "2" = "n = 40, p = 40",
                                                     "3" = "n = 10, p = 40")))) +
  labs(x = "Method", y = "Mean L1 accuracy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-1.png", sim.plot.1, width = 8, height = 5)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-1.png")

sim.plot.2 <- sim.res.df.2 %>% 
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
  labs(x = "Method", y = "M*") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-2.png", sim.plot.2, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-2.png")

sim.plot.3 <- sim.res.df.3 %>% 
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
  labs(x = "Method", y = "Mean lppd") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("Lasso/Lasso-plots/Lasso-sim-plot-3.png", sim.plot.3, width = 8, height = 3)
plot_crop("Lasso/Lasso-plots/Lasso-sim-plot-3.png")

sim.table.1 <- sim.res.df.2 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  mutate(n_log_cov_norm = -log(cov_norm + 0.00001)) %>% 
  dplyr::select(-c(mmd, cov_norm)) %>% 
  group_by(sim, method) %>% 
  summarise(mean_n_log_cov_norm = round(mean(n_log_cov_norm), 2)) %>% 
  pivot_wider(names_from = "sim", values_from = "mean_n_log_cov_norm")

sim.table.2 <- sim.res.df.4 %>% 
  mutate(method = factor(method, levels = c("mcmc", "mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  group_by(sim, method) %>% 
  summarise(mean_time = mean(time)) %>% 
  pivot_wider(names_from = "sim", values_from = "mean_time")

sim.table.3 <- sim.res.df.5 %>% 
  mutate(method = factor(method, levels = c("mcmc", "mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>%
  group_by(sim, iteration, method) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(sim, method) %>% 
  summarise(mean_max_r_hat = round(mean(max_r_hat), 2)) %>% 
  pivot_wider(names_from = "sim", values_from = "mean_max_r_hat")

## Benchmarks

bench.table.1 <- bench.res.df.1 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p"]])) ~ "beta",
                           TRUE ~ "kappa")) %>% 
  group_by(bench, block, method) %>% 
  summarise(mean_l1 = 100*round(mean(l1), 3)) %>% 
  pivot_wider(names_from = "bench", values_from = "mean_l1")

bench.table.2 <- bench.res.df.2 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  mutate(m_log_mmd = round(-log(mmd + 0.00001), 2)) %>% 
  dplyr::select(-c(mmd, cov_norm)) %>% 
  pivot_wider(names_from = "bench", values_from = "m_log_mmd")

bench.table.3 <- bench.res.df.3 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  filter(lppd != -Inf) %>% 
  group_by(bench, method) %>% 
  summarise(mean_lppd = round(mean(lppd), 2)) %>% 
  pivot_wider(names_from = "bench", values_from = "mean_lppd")

bench.table.4 <- bench.res.df.2 %>% 
  mutate(method = factor(method, levels = c("mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  mutate(m_log_cn = round(-log(cov_norm + 0.00001), 2)) %>% 
  dplyr::select(-c(mmd, cov_norm)) %>% 
  pivot_wider(names_from = "bench", values_from = "m_log_cn")

bench.table.5 <- bench.res.df.4 %>% 
  mutate(method = factor(method, levels = c("mcmc", "mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>% 
  group_by(bench, method) %>% 
  summarise(mean_time = round(mean(time), 2)) %>%
  pivot_wider(names_from = "bench", values_from = "mean_time")

bench.table.6 <- bench.res.df.5 %>% 
  mutate(method = factor(method, levels = c("mcmc", "mcmc-a", "mcmc-b", "mcmc-c", "ep", "mfvb"))) %>%
  group_by(bench, method) %>% 
  summarise(max_r_hat = round(max(r_hat), 2)) %>% 
  pivot_wider(names_from = "bench", values_from = "max_r_hat")
