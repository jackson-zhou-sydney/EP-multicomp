
# Plots and tables for lasso linear regression

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")
res.files <- list.files("Lasso/Results/")

## Simulations

### MCMC convergence tests

sim.conv.files <- res.files[grep("Simulations-conv-results", res.files)]
sim.r.hat.cdf <- data.frame()

for (file in sim.conv.files) {
  load(paste0("Lasso/Results/", file))
  sim.r.hat.cdf <- rbind(sim.r.hat.cdf, sim.r.hat.df)
}

sim.r.hat.table <- sim.r.hat.cdf %>% 
  group_by(sim, mcmc_iter, seed) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(sim, mcmc_iter) %>% 
  summarise(mean_max_r_hat = mean(max_r_hat))

### L1 accuracy

sim.l1.table <- sim.l1.df %>%
  mutate(block = case_when(j <= as.numeric(map(sim, ~sim.settings[[.]][["p"]])) ~ "beta",
                           TRUE ~ "kappa")) %>%
  group_by(seed, sim, iteration, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(seed, sim, method, block) %>%
  summarise(m_m_l1 = mean(m_l1)) %>%
  group_by(sim, method, block) %>%
  summarise(m_m_m_l1 = round(mean(m_m_l1), table.dp),
            sd_m_m_l1 = round(sd(m_m_l1), table.dp)) %>%
  arrange(sim, block, method)

### M-star

sim.m.star.table <- sim.mmd.df %>%
  mutate(m_star = -log(mmd + 0.00001)) %>%
  group_by(seed, sim, method) %>%
  summarise(m_m_star = mean(m_star)) %>%
  group_by(sim, method) %>%
  summarise(m_m_m_star = round(mean(m_m_star), table.dp),
            sd_m_m_star = round(sd(m_m_star), table.dp)) %>%
  arrange(sim, method)

### LPPD

sim.lppd.table <- sim.lppd.df %>%
  group_by(seed, sim, method) %>%
  summarise(m_lppd = mean(lppd)) %>%
  group_by(sim, method) %>%
  summarise(m_m_lppd = round(mean(m_lppd), table.dp),
            sd_m_lppd = round(sd(m_lppd), table.dp)) %>%
  arrange(sim, method)

### F-star

sim.f.star.table <- sim.cov.norm.df %>%
  mutate(f_star = -log(cov_norm + 0.00001)) %>%
  group_by(seed, sim, method) %>%
  summarise(m_f_star = mean(f_star)) %>%
  group_by(sim, method) %>%
  summarise(m_m_f_star = round(mean(m_f_star), table.dp),
            sd_m_f_star = round(sd(m_f_star), table.dp)) %>%
  arrange(sim, method)

### Run time

sim.time.table <- sim.time.df %>%
  group_by(seed, sim, method) %>%
  summarise(m_time = mean(time)) %>%
  group_by(sim, method) %>%
  summarise(m_m_time = round(mean(m_time), table.dp),
            sd_m_time = round(sd(m_time), table.dp)) %>%
  arrange(sim, method)

## Benchmarks

### MCMC convergence tests

bench.conv.files <- res.files[grep("Benchmarks-conv-results", res.files)]
bench.r.hat.cdf <- data.frame()

for (file in bench.conv.files) {
  load(paste0("Lasso/Results/", file))
  bench.r.hat.cdf <- rbind(bench.r.hat.cdf, bench.r.hat.df)
}

bench.r.hat.table <- bench.r.hat.cdf %>% 
  group_by(bench, mcmc_iter, seed) %>% 
  summarise(max_r_hat = max(r_hat)) %>% 
  group_by(bench, mcmc_iter) %>% 
  summarise(mean_max_r_hat = mean(max_r_hat))


