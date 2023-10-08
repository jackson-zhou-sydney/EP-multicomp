
# Big table

source("General-auxiliaries.R")

## Heteroscedastic linear regression

source("Hetero/Auxiliaries.R")
res.directory <- "Hetero/Results/"
plot.directory <- "Hetero/Plots/"
res.files <- list.files(res.directory)

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

bench.l1.table <- bench.l1.cdf %>%
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2"),
         l1 = 100*l1) %>% 
  group_by(seed, bench, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(bench, method, block) %>%
  summarise(m_m_l1 = round(mean(m_l1), l1.dp),
            sd_m_l1 = round(sd(m_l1), l1.dp)) %>%
  arrange(bench, block, method)

bench.m.star.table <- bench.mmd.cdf %>%
  mutate(m_star = -log(mmd + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_m_star = round(mean(m_star), table.dp),
            sd_m_star = round(sd(m_star), table.dp)) %>%
  arrange(bench, method)

bench.lppd.table <- bench.lppd.cdf %>%
  group_by(bench, method) %>%
  summarise(m_lppd = round(mean(lppd), table.dp),
            sd_lppd = round(sd(lppd), table.dp)) %>%
  arrange(bench, method)

bench.f.star.table <- bench.cov.norm.cdf %>%
  mutate(f_star = -log(cov_norm + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_f_star = round(mean(f_star), table.dp),
            sd_f_star = round(sd(f_star), table.dp)) %>%
  arrange(bench, method)

bench.time.table <- bench.time.cdf %>%
  group_by(bench, method) %>%
  summarise(m_time = round(mean(time), table.dp),
            sd_time = round(sd(time), table.dp)) %>%
  arrange(bench, method)

bench.l1.table.clean <- bench.l1.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.1f", m_m_l1), " $\\pm$ ", sprintf("%3.1f", sd_m_l1))) %>% 
  select(-m_m_l1, -sd_m_l1) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  ungroup()

hetero.bench.beta.1 <- bench.l1.table.clean %>% 
  filter(block == "beta_1") %>% 
  select(-block) %>% 
  add_row(method = factor("MFVB"), `4` = "---") %>% 
  arrange(method)

hetero.bench.beta.2 <- bench.l1.table.clean %>% 
  filter(block == "beta_2") %>% 
  select(-block) %>% 
  add_row(method = factor("MFVB"), `4` = "---") %>% 
  arrange(method)

hetero.bench.m.star <- bench.m.star.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.2f", m_m_star), " $\\pm$ ", sprintf("%3.2f", sd_m_star))) %>% 
  select(-m_m_star, -sd_m_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("MFVB"), `4` = "---") %>% 
  arrange(method)

hetero.bench.lppd <- bench.lppd.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%4.2f", m_lppd), " $\\pm$ ", sprintf("%3.2f", sd_lppd))) %>% 
  select(-m_lppd, -sd_lppd) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("MFVB"), `4` = "---") %>% 
  arrange(method)

hetero.bench.f.star <- bench.f.star.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.2f", m_f_star), " $\\pm$ ", sprintf("%3.2f", sd_f_star))) %>% 
  select(-m_f_star, -sd_f_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("MFVB"), `4` = "---") %>% 
  arrange(method)

hetero.bench.time <- bench.time.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%4.2f", m_time), " $\\pm$ ", sprintf("%3.2f", sd_time))) %>% 
  select(-m_time, -sd_time) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("MFVB"), `4` = "---") %>% 
  arrange(method)

## Lasso linear regression

source("Lasso/Auxiliaries.R")
res.directory <- "Lasso/Results/"
plot.directory <- "Lasso/Plots/"
res.files <- list.files(res.directory)

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

bench.l1.table <- bench.l1.cdf %>%
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p"]])) ~ "beta",
                           TRUE ~ "kappa"),
         l1 = 100*l1) %>%
  group_by(seed, bench, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(bench, method, block) %>%
  summarise(m_m_l1 = round(mean(m_l1), l1.dp),
            sd_m_l1 = round(sd(m_l1), l1.dp)) %>%
  arrange(bench, block, method)

bench.m.star.table <- bench.mmd.cdf %>%
  mutate(m_star = -log(mmd + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_m_star = round(mean(m_star), table.dp),
            sd_m_star = round(sd(m_star), table.dp)) %>%
  arrange(bench, method)

bench.lppd.table <- bench.lppd.cdf %>%
  group_by(bench, method) %>%
  summarise(m_lppd = round(mean(lppd), table.dp),
            sd_lppd = round(sd(lppd), table.dp)) %>%
  arrange(bench, method)

bench.f.star.table <- bench.cov.norm.cdf %>%
  mutate(f_star = -log(cov_norm + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_f_star = round(mean(f_star), table.dp),
            sd_f_star = round(sd(f_star), table.dp)) %>%
  arrange(bench, method)

bench.time.table <- bench.time.cdf %>%
  group_by(bench, method) %>%
  summarise(m_time = round(mean(time), table.dp),
            sd_time = round(sd(time), table.dp)) %>%
  arrange(bench, method)

bench.l1.table.clean <- bench.l1.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.1f", m_m_l1), " $\\pm$ ", sprintf("%3.1f", sd_m_l1))) %>% 
  select(-m_m_l1, -sd_m_l1) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  ungroup()

lasso.bench.beta <- bench.l1.table.clean %>% 
  filter(block == "beta") %>% 
  select(-block) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

lasso.bench.kappa <- bench.l1.table.clean %>% 
  filter(block == "kappa") %>% 
  select(-block) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

lasso.bench.m.star <- bench.m.star.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.2f", m_m_star), " $\\pm$ ", sprintf("%3.2f", sd_m_star))) %>% 
  select(-m_m_star, -sd_m_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

lasso.bench.lppd <- bench.lppd.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%4.2f", m_lppd), " $\\pm$ ", sprintf("%3.2f", sd_lppd))) %>% 
  select(-m_lppd, -sd_lppd) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

lasso.bench.f.star <- bench.f.star.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.2f", m_f_star), " $\\pm$ ", sprintf("%3.2f", sd_f_star))) %>% 
  select(-m_f_star, -sd_f_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

lasso.bench.time <- bench.time.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%4.2f", m_time), " $\\pm$ ", sprintf("%3.2f", sd_time))) %>% 
  select(-m_time, -sd_time) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

## Quantile regression

source("Quantile/Auxiliaries.R")
res.directory <- "Quantile/Results/"
plot.directory <- "Quantile/Plots/"
res.files <- list.files(res.directory)

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

bench.l1.table <- bench.l1.cdf %>%
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p"]])) ~ "beta",
                           TRUE ~ "kappa"),
         l1 = 100*l1) %>%
  group_by(seed, bench, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(bench, method, block) %>%
  summarise(m_m_l1 = round(mean(m_l1), l1.dp),
            sd_m_l1 = round(sd(m_l1), l1.dp)) %>%
  arrange(bench, block, method)

bench.m.star.table <- bench.mmd.cdf %>%
  mutate(m_star = -log(mmd + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_m_star = round(mean(m_star), table.dp),
            sd_m_star = round(sd(m_star), table.dp)) %>%
  arrange(bench, method)

bench.lppd.table <- bench.lppd.cdf %>%
  group_by(bench, method) %>%
  summarise(m_lppd = round(mean(lppd), table.dp),
            sd_lppd = round(sd(lppd), table.dp)) %>%
  arrange(bench, method)

bench.f.star.table <- bench.cov.norm.cdf %>%
  mutate(f_star = -log(cov_norm + log.lb)) %>%
  group_by(bench, method) %>%
  summarise(m_f_star = round(mean(f_star), table.dp),
            sd_f_star = round(sd(f_star), table.dp)) %>%
  arrange(bench, method)

bench.time.table <- bench.time.cdf %>%
  group_by(bench, method) %>%
  summarise(m_time = round(mean(time), table.dp),
            sd_time = round(sd(time), table.dp)) %>%
  arrange(bench, method)

bench.l1.table.clean <- bench.l1.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.1f", m_m_l1), " $\\pm$ ", sprintf("%3.1f", sd_m_l1))) %>% 
  select(-m_m_l1, -sd_m_l1) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  ungroup()

quantile.bench.beta <- bench.l1.table.clean %>% 
  filter(block == "beta") %>% 
  select(-block) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

quantile.bench.kappa <- bench.l1.table.clean %>% 
  filter(block == "kappa") %>% 
  select(-block) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

quantile.bench.m.star <- bench.m.star.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.2f", m_m_star), " $\\pm$ ", sprintf("%3.2f", sd_m_star))) %>% 
  select(-m_m_star, -sd_m_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

quantile.bench.lppd <- bench.lppd.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%4.2f", m_lppd), " $\\pm$ ", sprintf("%3.2f", sd_lppd))) %>% 
  select(-m_lppd, -sd_lppd) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

quantile.bench.f.star <- bench.f.star.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%3.2f", m_f_star), " $\\pm$ ", sprintf("%3.2f", sd_f_star))) %>% 
  select(-m_f_star, -sd_f_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

quantile.bench.time <- bench.time.table %>% 
  filter(bench == 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "mfvb" ~ "MFVB"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C", "MFVB")),
         res = paste0(sprintf("%4.2f", m_time), " $\\pm$ ", sprintf("%3.2f", sd_time))) %>% 
  select(-m_time, -sd_time) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  add_row(method = factor("Laplace"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-A"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-B"), `4` = "---") %>% 
  add_row(method = factor("Pathfinder-C"), `4` = "---") %>% 
  arrange(method)

## Combining

chunk.1 <- cbind(hetero.bench.beta.1, lasso.bench.beta, quantile.bench.beta,
      hetero.bench.beta.2, lasso.bench.kappa, quantile.bench.kappa) %>% 
  `colnames<-`(letters[1:12]) %>% 
  select(-c("c", "e", "g", "i", "k")) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

chunk.2 <- cbind(hetero.bench.m.star, lasso.bench.m.star, quantile.bench.m.star,
                 hetero.bench.lppd, lasso.bench.lppd, quantile.bench.lppd) %>% 
  `colnames<-`(letters[1:12]) %>% 
  select(-c("c", "e", "g", "i", "k")) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

chunk.3 <- cbind(hetero.bench.f.star, lasso.bench.f.star, quantile.bench.f.star,
                 hetero.bench.time, lasso.bench.time, quantile.bench.time) %>% 
  `colnames<-`(letters[1:12]) %>% 
  select(-c("c", "e", "g", "i", "k")) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

big.latex <- paste0("\\begin{table}\n",
                    "\\centering\n",
                    "\\resizebox{\\linewidth}{!}{\n",
                    "\\begin{tabular}{@{}ccccccc@{}}\n", 
                    "\\toprule\n",
                    "& Hetero. & Lasso & Quantile & Hetero. & Lasso & Quantile \\\\ \\midrule\n",
                    "& $\\vbeta_1$ $L^1$ acc. & $\\vbeta$ $L^1$ acc. & $\\vbeta$ $L^1$ acc. & $\\vbeta_2$ $L^1$ acc. & $\\kappa$ $L^1$ acc. & $\\kappa$ $L^1$ acc. \\\\\n",
                    chunk.1, "\\\\ \\midrule\n",
                    "& \\multicolumn{3}{c}{$M^*$} & \\multicolumn{3}{c}{lppd} \\\\\n",
                    chunk.2, "\\\\ \\midrule\n",
                    "& \\multicolumn{3}{c}{$F^*$} & \\multicolumn{3}{c}{Run time (seconds)} \\\\\n",
                    chunk.3, "\\\\ \\bottomrule\n",
                    "\\end{tabular}}\n",
                    "\\caption{Results (mean $\\pm$ SD) for energy dataset. Each cell represents eight repetitions.}\n",
                    "\\label{table:big}\n",
                    "\\end{table}")

cat(big.latex)
