
# Plots and tables for heteroscedastic linear regression

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")
res.directory <- "Hetero/Results/"
res.files <- list.files(res.directory)

## Simulations

sim.conv.files <- res.files[grep("Simulations-conv-results", res.files)]
sim.r.hat.cdf <- data.frame()

for (file in sim.conv.files) {
  load(paste0(res.directory, file))
  sim.r.hat.cdf <- rbind(sim.r.hat.cdf, sim.r.hat.df)
}

sim.files <- res.files[grep("Simulations-results", res.files)]
sim.files <- sim.files[-grep("MCMC-G", sim.files)]
sim.l1.cdf <- data.frame()
sim.mmd.cdf <- data.frame()
sim.lppd.cdf <- data.frame()
sim.cov.norm.cdf <- data.frame()
sim.time.cdf <- data.frame()

for (file in sim.files) {
  load(paste0(res.directory, file))
  sim.l1.cdf <- rbind(sim.l1.cdf, sim.l1.df)
  sim.mmd.cdf <- rbind(sim.mmd.cdf, sim.mmd.df)
  sim.lppd.cdf <- rbind(sim.lppd.cdf, sim.lppd.df)
  sim.cov.norm.cdf <- rbind(sim.cov.norm.cdf, sim.cov.norm.df)
  sim.time.cdf <- rbind(sim.time.cdf, sim.time.df)
}

### MCMC convergence tests

sim.r.hat.table <- sim.r.hat.cdf %>% 
  group_by(sim, mcmc_iter, seed) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(sim, mcmc_iter) %>% 
  summarise(mean_max_r_hat = mean(max_r_hat))

save(sim.r.hat.table, file = paste0(res.directory, "Simulations-conv-table.RData"))

### L1 accuracy

sim.l1.plot <- sim.l1.cdf %>% 
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  mutate(block = case_when(j <= as.numeric(map(sim, ~sim.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2")) %>% 
  group_by(seed, sim, iteration, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(seed, sim, method, block) %>%
  summarise(m_m_l1 = mean(m_l1)) %>% 
  ggplot(mapping = aes(x = method, y = m_m_l1)) +
  geom_boxplot() +
  ggh4x::facet_grid2(block ~ sim, scales = "free_y", independent = "y",
                     labeller = labeller(block = as_labeller(c("beta_1" = "Beta_1", "beta_2" = "Beta_2")),
                                         sim = as_labeller(sim.labels))) +
  labs(x = "Method", y = "Mean L1 accuracy across iterations and marginals") +
  theme_bw()

sim.l1.table <- sim.l1.cdf %>%
  mutate(block = case_when(j <= as.numeric(map(sim, ~sim.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2")) %>% 
  group_by(seed, sim, iteration, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(seed, sim, method, block) %>%
  summarise(m_m_l1 = mean(m_l1)) %>%
  group_by(sim, method, block) %>%
  summarise(m_m_m_l1 = round(mean(m_m_l1), table.dp),
            sd_m_m_l1 = round(sd(m_m_l1), table.dp)) %>%
  arrange(sim, block, method)

### M-star

sim.m.star.plot <- sim.mmd.cdf %>% 
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  mutate(m_star = -log(mmd + log.lb)) %>%
  group_by(seed, sim, method) %>%
  summarise(m_m_star = mean(m_star)) %>% 
  ggplot(mapping = aes(x = method, y = m_m_star)) +
  geom_boxplot() +
  facet_wrap(~sim, scales = "free_y",
             labeller = as_labeller(sim.labels)) +
  labs(x = "Method", y = "Mean M-star across iterations") +
  theme_bw()

sim.m.star.table <- sim.mmd.cdf %>%
  mutate(m_star = -log(mmd + log.lb)) %>%
  group_by(seed, sim, method) %>%
  summarise(m_m_star = mean(m_star)) %>%
  group_by(sim, method) %>%
  summarise(m_m_m_star = round(mean(m_m_star), table.dp),
            sd_m_m_star = round(sd(m_m_star), table.dp)) %>%
  arrange(sim, method)

### LPPD

sim.lppd.plot <- sim.lppd.cdf %>% 
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  group_by(seed, sim, method) %>%
  summarise(m_lppd = mean(lppd)) %>% 
  ggplot(mapping = aes(x = method, y = m_lppd)) +
  geom_boxplot() +
  facet_wrap(~sim, scales = "free_y",
             labeller = as_labeller(sim.labels)) +
  labs(x = "Method", y = "Mean lppd across iterations") +
  theme_bw()

sim.lppd.table <- sim.lppd.cdf %>%
  group_by(seed, sim, method) %>%
  summarise(m_lppd = mean(lppd)) %>%
  group_by(sim, method) %>%
  summarise(m_m_lppd = round(mean(m_lppd), table.dp),
            sd_m_lppd = round(sd(m_lppd), table.dp)) %>%
  arrange(sim, method)

### F-star

sim.f.star.plot <- sim.cov.norm.cdf %>% 
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  mutate(f_star = -log(cov_norm + log.lb)) %>%
  group_by(seed, sim, method) %>%
  summarise(m_f_star = mean(f_star)) %>% 
  ggplot(mapping = aes(x = method, y = m_f_star)) +
  geom_boxplot() +
  facet_wrap(~sim, scales = "free_y",
             labeller = as_labeller(sim.labels)) +
  labs(x = "Method", y = "Mean F-star across iterations") +
  theme_bw()

sim.f.star.table <- sim.cov.norm.cdf %>%
  mutate(f_star = -log(cov_norm + log.lb)) %>%
  group_by(seed, sim, method) %>%
  summarise(m_f_star = mean(f_star)) %>%
  group_by(sim, method) %>%
  summarise(m_m_f_star = round(mean(m_f_star), table.dp),
            sd_m_f_star = round(sd(m_f_star), table.dp)) %>%
  arrange(sim, method)

### Run time

sim.time.plot <- sim.time.cdf %>% 
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  group_by(seed, sim, method) %>%
  summarise(m_time = mean(time)) %>% 
  ggplot(mapping = aes(x = method, y = m_time)) +
  geom_boxplot() +
  facet_wrap(~sim, scales = "free_y",
             labeller = as_labeller(sim.labels)) +
  scale_y_log10() +
  labs(x = "Method", y = "Mean time across iterations (log scale)") +
  theme_bw()

sim.time.table <- sim.time.cdf %>%
  group_by(seed, sim, method) %>%
  summarise(m_time = mean(time)) %>%
  group_by(sim, method) %>%
  summarise(m_m_time = round(mean(m_time), table.dp),
            sd_m_time = round(sd(m_time), table.dp)) %>%
  arrange(sim, method)

## Benchmarks

bench.conv.files <- res.files[grep("Benchmarks-conv-results", res.files)]
bench.r.hat.cdf <- data.frame()

for (file in bench.conv.files) {
  load(paste0(res.directory, file))
  bench.r.hat.cdf <- rbind(bench.r.hat.cdf, bench.r.hat.df)
}

bench.files <- res.files[grep("Benchmarks-results", res.files)]
bench.files <- bench.files[-grep("MCMC-G", bench.files)]
big.files <- res.files[grep("Big-results", res.files)]
bench.l1.cdf <- data.frame()
bench.mmd.cdf <- data.frame()
bench.lppd.cdf <- data.frame()
bench.cov.norm.cdf <- data.frame()
bench.time.cdf <- data.frame()

for (file in c(bench.files, big.files)) {
  load(paste0(res.directory, file))
  bench.l1.cdf <- rbind(bench.l1.cdf, bench.l1.df)
  bench.mmd.cdf <- rbind(bench.mmd.cdf, bench.mmd.df)
  bench.lppd.cdf <- rbind(bench.lppd.cdf, bench.lppd.df)
  bench.cov.norm.cdf <- rbind(bench.cov.norm.cdf, bench.cov.norm.df)
  bench.time.cdf <- rbind(bench.time.cdf, bench.time.df)
}

### MCMC convergence tests

bench.r.hat.table <- bench.r.hat.cdf %>% 
  group_by(bench, mcmc_iter, seed) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(bench, mcmc_iter) %>% 
  summarise(mean_max_r_hat = mean(max_r_hat))

save(bench.r.hat.table, file = paste0(res.directory, "Benchmarks-conv-table.RData"))

### L1 accuracy

bench.l1.plot <- bench.l1.cdf %>%
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2")) %>% 
  group_by(seed, bench, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  ggplot(mapping = aes(x = method, y = m_l1)) +
  geom_boxplot() +
  ggh4x::facet_grid2(block ~ bench, scales = "free_y", independent = "y",
                     labeller = labeller(block = as_labeller(c("beta_1" = "Beta_1", "beta_2" = "Beta_2")),
                                         bench = as_labeller(bench.labels))) +
  labs(x = "Method", y = "Mean L1 accuracy across marginals") +
  theme_bw()

bench.l1.table <- bench.l1.cdf %>%
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2")) %>% 
  group_by(seed, bench, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(bench, method, block) %>%
  summarise(m_m_l1 = round(mean(m_l1), table.dp),
            sd_m_l1 = round(sd(m_l1), table.dp)) %>%
  arrange(bench, block, method)

### M-star
  
bench.m.star.plot <- bench.mmd.cdf %>%
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  mutate(m_star = -log(mmd + log.lb)) %>% 
  ggplot(mapping = aes(x = method, y = m_star)) +
  geom_boxplot() +
  facet_wrap(~bench, scales = "free_y",
             labeller = as_labeller(bench.labels)) +
  labs(x = "Method", y = "M-star") +
  theme_bw()

bench.m.star.table <- bench.mmd.cdf %>%
  mutate(m_star = -log(mmd + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_m_star = round(mean(m_star), table.dp),
            sd_m_star = round(sd(m_star), table.dp)) %>%
  arrange(bench, method)

### LPPD

bench.lppd.plot <- bench.lppd.cdf %>% 
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  ggplot(mapping = aes(x = method, y = lppd)) +
  geom_boxplot() +
  facet_wrap(~bench, scales = "free_y",
             labeller = as_labeller(bench.labels)) +
  labs(x = "Method", y = "lppd") +
  theme_bw()

bench.lppd.table <- bench.lppd.cdf %>%
  group_by(bench, method) %>%
  summarise(m_lppd = round(mean(lppd), table.dp),
            sd_lppd = round(sd(lppd), table.dp)) %>%
  arrange(bench, method)

### F-star

bench.f.star.plot <- bench.cov.norm.cdf %>%
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  mutate(f_star = -log(cov_norm + log.lb)) %>% 
  ggplot(mapping = aes(x = method, y = f_star)) +
  geom_boxplot() +
  facet_wrap(~bench, scales = "free_y",
             labeller = as_labeller(bench.labels)) +
  labs(x = "Method", y = "F-star") +
  theme_bw()

bench.f.star.table <- bench.cov.norm.cdf %>%
  mutate(f_star = -log(cov_norm + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_f_star = round(mean(f_star), table.dp),
            sd_f_star = round(sd(f_star), table.dp)) %>%
  arrange(bench, method)

### Run time

bench.time.plot <- bench.time.cdf %>% 
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  ggplot(mapping = aes(x = method, y = time)) +
  geom_boxplot() +
  facet_wrap(~bench, scales = "free_y",
             labeller = as_labeller(bench.labels)) +
  scale_y_log10() +
  labs(x = "Method", y = "Time (log scale)") +
  theme_bw()

bench.time.table <- bench.time.cdf %>%
  group_by(bench, method) %>%
  summarise(m_time = round(mean(time), table.dp),
            sd_time = round(sd(time), table.dp)) %>%
  arrange(bench, method)
