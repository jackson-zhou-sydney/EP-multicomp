
# Plots and tables for heteroscedastic linear regression

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")
res.directory <- "Hetero/Results/"
plot.directory <- "Hetero/Plots/"
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
                           TRUE ~ "beta_2"),
         l1 = 100*l1) %>% 
  group_by(seed, sim, iteration, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(seed, sim, method, block) %>%
  summarise(m_m_l1 = mean(m_l1)) %>% 
  ggplot(mapping = aes(x = method, y = m_m_l1)) +
  geom_boxplot() +
  ggh4x::facet_grid2(block ~ sim, scales = "free_y", independent = "y",
                     labeller = labeller(sim = as_labeller(sim.labels))) +
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  labs(x = "Method", y = "Mean L1 accuracy across iterations and marginals") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

sim.l1.table <- sim.l1.cdf %>%
  mutate(block = case_when(j <= as.numeric(map(sim, ~sim.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2"),
         l1 = 100*l1) %>% 
  group_by(seed, sim, iteration, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  group_by(seed, sim, method, block) %>%
  summarise(m_m_l1 = mean(m_l1)) %>%
  group_by(sim, method, block) %>%
  summarise(m_m_m_l1 = round(mean(m_m_l1), l1.dp),
            sd_m_m_l1 = round(sd(m_m_l1), l1.dp)) %>%
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
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  labs(x = "Method", y = "Mean M-star across iterations") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  labs(x = "Method", y = "Mean lppd across iterations") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  labs(x = "Method", y = "Mean F-star across iterations") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  scale_y_log10() +
  labs(x = "Method", y = "Mean time across iterations (log scale)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

sim.time.table <- sim.time.cdf %>%
  group_by(seed, sim, method) %>%
  summarise(m_time = mean(time)) %>%
  group_by(sim, method) %>%
  summarise(m_m_time = round(mean(m_time), table.dp),
            sd_m_time = round(sd(m_time), table.dp)) %>%
  arrange(sim, method)

### Combined LaTeX table

sim.l1.table.clean <- sim.l1.table %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%3.1f", m_m_m_l1), " $\\pm$ ", sprintf("%3.1f", sd_m_m_l1))) %>% 
  select(-m_m_m_l1, -sd_m_m_l1) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = sim, values_from = res) %>% 
  ungroup()

sim.beta.1 <- sim.l1.table.clean %>% 
  filter(block == "beta_1") %>% 
  select(-block)

sim.beta.2 <- sim.l1.table.clean %>% 
  filter(block == "beta_2") %>% 
  select(-block)

sim.m.star <- sim.m.star.table %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%3.2f", m_m_m_star), " $\\pm$ ", sprintf("%3.2f", sd_m_m_star))) %>% 
  select(-m_m_m_star, -sd_m_m_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = sim, values_from = res)

sim.lppd <- sim.lppd.table %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%4.2f", m_m_lppd), " $\\pm$ ", sprintf("%3.2f", sd_m_lppd))) %>% 
  select(-m_m_lppd, -sd_m_lppd) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = sim, values_from = res)

sim.f.star <- sim.f.star.table %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%3.2f", m_m_f_star), " $\\pm$ ", sprintf("%3.2f", sd_m_f_star))) %>% 
  select(-m_m_f_star, -sd_m_f_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = sim, values_from = res)

sim.time <- sim.time.table %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%4.2f", m_m_time), " $\\pm$ ", sprintf("%3.2f", sd_m_time))) %>% 
  select(-m_m_time, -sd_m_time) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = sim, values_from = res)

chunk.1 <- cbind(sim.beta.1, sim.beta.2) %>% 
  `colnames<-`(letters[1:8]) %>% 
  select(-e) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

chunk.2 <- cbind(sim.m.star, sim.lppd) %>% 
  `colnames<-`(letters[1:8]) %>% 
  select(-e) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

chunk.3 <- cbind(sim.f.star, sim.time) %>% 
  `colnames<-`(letters[1:8]) %>% 
  select(-e) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

sim.latex <- paste0("\\begin{table}\n",
                    "\\centering\n",
                    "\\resizebox{\\linewidth}{!}{\n",
                    "\\begin{tabular}{@{}ccccccc@{}}\n", 
                    "\\toprule\n",
                    "& Setting 1 & Setting 2 & Setting 3 & Setting 1 & Setting 2 & Setting 3 \\\\ \\midrule\n",
                    "& \\multicolumn{3}{c}{$\\vbeta_1$ $L^1$ accuracy} & \\multicolumn{3}{c}{$\\vbeta_2$ $L^1$ accuracy} \\\\\n",
                    chunk.1, "\\\\ \\midrule\n",
                    "& \\multicolumn{3}{c}{$M^*$} & \\multicolumn{3}{c}{lppd} \\\\\n",
                    chunk.2, "\\\\ \\midrule\n",
                    "& \\multicolumn{3}{c}{$F^*$} & \\multicolumn{3}{c}{Run time (seconds)} \\\\\n",
                    chunk.3, "\\\\ \\bottomrule\n",
                    "\\end{tabular}}\n",
                    "\\caption{Results (mean $\\pm$ SD) across heteroscedastic linear regression simulations. Each cell represents thirty repetitions. Settings are organized by order of appearance in the main text.}\n",
                    "\\label{table:heterosim}\n",
                    "\\end{table}")

cat(sim.latex)

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
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "LM", "GVB-A", "GVB-B", "GVB-C"))) %>% 
  mutate(block = case_when(j <= as.numeric(map(bench, ~bench.settings[[.]][["p.1"]])) ~ "beta_1",
                           TRUE ~ "beta_2"),
         l1 = 100*l1) %>% 
  group_by(seed, bench, method, block) %>% 
  summarise(m_l1 = mean(l1)) %>%
  ggplot(mapping = aes(x = method, y = m_l1)) +
  geom_boxplot() +
  ggh4x::facet_grid2(block ~ bench, scales = "free_y", independent = "y",
                     labeller = labeller(bench = as_labeller(c("1" = "Food",
                                                               "2" = "Salary",
                                                               "3" = "Sniffer",
                                                               "4" = "Energy")))) +
  scale_x_discrete(labels = c("MCMC" = "ML", "MCMC-S" = "MS", "EP" = "EP", "EP-2D" = "E2", "GVB-A" = "PA", "GVB-B" = "PB", "GVB-C" = "PC", "LM" = "LA")) +
  labs(x = "Method", y = "Mean L1 accuracy") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(paste0(plot.directory, "Benchmarks-L1.png"), plot = bench.l1.plot, dpi = 600, width = 15, height = 8, units = "cm")

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

### M-star
  
bench.m.star.plot <- bench.mmd.cdf %>%
  mutate(method = factor(toupper(method), levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "GVB-A", "GVB-B", "GVB-C", "LM"))) %>% 
  mutate(m_star = -log(mmd + log.lb)) %>% 
  ggplot(mapping = aes(x = method, y = m_star)) +
  geom_boxplot() +
  facet_wrap(~bench, scales = "free_y",
             labeller = as_labeller(bench.labels)) +
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  labs(x = "Method", y = "M-star") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  labs(x = "Method", y = "lppd") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  labs(x = "Method", y = "F-star") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

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
  scale_x_discrete(labels = c("MCMC" = "MCMC", "MCMC-S" = "MCMC-S", "EP" = "EP", "EP-2D" = "EP-2D", "GVB-A" = "Pathfinder-A", "GVB-B" = "Pathfinder-B", "GVB-C" = "Pathfinder-C", "LM" = "Laplace")) +
  scale_y_log10() +
  labs(x = "Method", y = "Time (log scale)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

bench.time.table <- bench.time.cdf %>%
  group_by(bench, method) %>%
  summarise(m_time = round(mean(time), table.dp),
            sd_time = round(sd(time), table.dp)) %>%
  arrange(bench, method)

### Accuracy vs. time

mean.l1.beta.1 <- bench.l1.table %>% 
  ungroup() %>% 
  filter(bench == 4) %>% 
  filter(block == "beta_1") %>% 
  select(method, m_m_l1)

mean.l1.beta.2 <- bench.l1.table %>% 
  ungroup() %>% 
  filter(bench == 4) %>% 
  filter(block == "beta_2") %>% 
  select(method, m_m_l1)

mean.time <- bench.time.table %>% 
  ungroup() %>% 
  filter(bench == 4) %>% 
  select(method, m_time)

beta.1.plot <- merge(mean.l1.beta.1, mean.time, by = "method") %>% 
  mutate(method_clean = case_when(method == "mcmc" ~ "ML",
                                  method == "mcmc-s" ~ "MS",
                                  method == "ep" ~ "EP",
                                  method == "ep-2d" ~ "E2",
                                  method == "gvb-a" ~ "PA",
                                  method == "gvb-b" ~ "PB",
                                  method == "gvb-c" ~ "PC",
                                  method == "lm" ~ "LA")) %>% 
  ggplot(aes(x = m_time, y = m_m_l1, label = method_clean)) +
  geom_point() +
  geom_text(size = 4, vjust = 0.5, hjust = -0.3) +
  lims(y = c(-4, 100)) +
  scale_x_log10() +
  coord_cartesian(xlim = c(1, 600000)) +
  labs(x = "Run time (seconds)",
       y = "beta_1 mean L1 accuracy") +
  theme_bw()

ggsave(paste0(plot.directory, "Benchmarks-beta-1.png"), plot = beta.1.plot, dpi = 600, width = 15, height = 7.5, units = "cm")

beta.2.plot <- merge(mean.l1.beta.2, mean.time, by = "method") %>% 
  mutate(method_clean = case_when(method == "mcmc" ~ "ML",
                                  method == "mcmc-s" ~ "MS",
                                  method == "ep" ~ "EP",
                                  method == "ep-2d" ~ "E2",
                                  method == "gvb-a" ~ "PA",
                                  method == "gvb-b" ~ "PB",
                                  method == "gvb-c" ~ "PC",
                                  method == "lm" ~ "LA"),
         vjust = case_when(method == "gvb-b" ~ 0.75,
                           method == "gvb-c" ~ 0.25,
                           TRUE ~ 0.5)) %>% 
  ggplot(aes(x = m_time, y = m_m_l1, label = method_clean)) +
  geom_point() +
  geom_text(aes(vjust = vjust), size = 4, hjust = -0.3) +
  lims(y = c(-4, 100)) +
  scale_x_log10() +
  coord_cartesian(xlim = c(1, 600000)) +
  labs(x = "Run time (seconds)",
       y = "beta_2 mean L1 accuracy") +
  theme_bw()

ggsave(paste0(plot.directory, "Benchmarks-beta-2.png"), plot = beta.2.plot, dpi = 600, width = 15, height = 7.5, units = "cm")

### Combined LaTeX table

bench.l1.table.clean <- bench.l1.table %>% 
  filter(bench != 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%3.1f", m_m_l1), " $\\pm$ ", sprintf("%3.1f", sd_m_l1))) %>% 
  select(-m_m_l1, -sd_m_l1) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res) %>% 
  ungroup()

bench.beta.1 <- bench.l1.table.clean %>% 
  filter(block == "beta_1") %>% 
  select(-block)

bench.beta.2 <- bench.l1.table.clean %>% 
  filter(block == "beta_2") %>% 
  select(-block)

bench.m.star <- bench.m.star.table %>% 
  filter(bench != 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%3.2f", m_m_star), " $\\pm$ ", sprintf("%3.2f", sd_m_star))) %>% 
  select(-m_m_star, -sd_m_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res)

bench.lppd <- bench.lppd.table %>% 
  filter(bench != 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%4.2f", m_lppd), " $\\pm$ ", sprintf("%3.2f", sd_lppd))) %>% 
  select(-m_lppd, -sd_lppd) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res)

bench.f.star <- bench.f.star.table %>% 
  filter(bench != 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%3.2f", m_f_star), " $\\pm$ ", sprintf("%3.2f", sd_f_star))) %>% 
  select(-m_f_star, -sd_f_star) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res)

bench.time <- bench.time.table %>% 
  filter(bench != 4) %>% 
  mutate(method = case_when(method == "mcmc" ~ "MCMC",
                            method == "mcmc-s" ~ "MCMC-S",
                            method == "ep" ~ "EP",
                            method == "ep-2d" ~ "EP-2D",
                            method == "gvb-a" ~ "Pathfinder-A",
                            method == "gvb-b" ~ "Pathfinder-B",
                            method == "gvb-c" ~ "Pathfinder-C",
                            method == "lm" ~ "Laplace"),
         method = factor(method, levels = c("MCMC", "MCMC-S", "EP", "EP-2D", "Laplace", "Pathfinder-A", "Pathfinder-B", "Pathfinder-C")),
         res = paste0(sprintf("%4.2f", m_time), " $\\pm$ ", sprintf("%3.2f", sd_time))) %>% 
  select(-m_time, -sd_time) %>% 
  arrange(method) %>% 
  pivot_wider(names_from = bench, values_from = res)

chunk.1 <- cbind(bench.beta.1, bench.beta.2) %>% 
  `colnames<-`(letters[1:8]) %>% 
  select(-e) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

chunk.2 <- cbind(bench.m.star, bench.lppd) %>% 
  `colnames<-`(letters[1:8]) %>% 
  select(-e) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

chunk.3 <- cbind(bench.f.star, bench.time) %>% 
  `colnames<-`(letters[1:8]) %>% 
  select(-e) %>% 
  mutate(a = as.character(a)) %>% 
  apply(1, function(x) paste0(x, collapse = " & ")) %>% 
  paste0(collapse = " \\\\\n")

bench.latex <- paste0("\\begin{table}\n",
                      "\\centering\n",
                      "\\resizebox{\\linewidth}{!}{\n",
                      "\\begin{tabular}{@{}ccccccc@{}}\n", 
                      "\\toprule\n",
                      "& Food & Salary & Sniffer & Food & Salary & Sniffer \\\\ \\midrule\n",
                      "& \\multicolumn{3}{c}{$\\vbeta_1$ $L^1$ accuracy} & \\multicolumn{3}{c}{$\\vbeta_2$ $L^1$ accuracy} \\\\\n",
                      chunk.1, "\\\\ \\midrule\n",
                      "& \\multicolumn{3}{c}{$M^*$} & \\multicolumn{3}{c}{lppd} \\\\\n",
                      chunk.2, "\\\\ \\midrule\n",
                      "& \\multicolumn{3}{c}{$F^*$} & \\multicolumn{3}{c}{Run time (seconds)} \\\\\n",
                      chunk.3, "\\\\ \\bottomrule\n",
                      "\\end{tabular}}\n",
                      "\\caption{Results (mean $\\pm$ SD) across heteroscedastic linear regression benchmarks. Each cell represents thirty repetitions.}\n",
                      "\\label{table:heterobench}\n",
                      "\\end{table}")

cat(bench.latex)
