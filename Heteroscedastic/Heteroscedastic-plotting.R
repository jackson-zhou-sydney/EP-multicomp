
# Plots for heteroscedastic linear regression

source("EP-general-auxiliaries.R")
source("Heteroscedastic/Heteroscedastic-auxiliaries.R")

## Simulations

sim.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Heteroscedastic/Heteroscedastic-data/Sim-", type.iter, "-iter-01.RData"))
  p.1 <- ncol(X.1)
  
  load(paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-Laplace.RData"))
  laplace.df <- results.df %>% mutate(method = "laplace")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Sim-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, laplace.df, ep.df) %>% 
    filter(j <= p.1) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_laplace = (mu_laplace - mu_mcmc)/sigma_2_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/sigma_2_mcmc,
           sigma2_laplace = (sigma_2_laplace - sigma_2_mcmc)/sigma_2_mcmc,
           sigma2_ep = (sigma_2_ep - sigma_2_mcmc)/sigma_2_mcmc) %>% 
    select(iteration, j, mu_laplace, mu_ep, sigma2_laplace, sigma2_ep) %>% 
    pivot_longer(cols = c(mu_laplace, mu_ep, sigma2_laplace, sigma2_ep)) %>% 
    separate(col = name, into = c("statistic", "method")) %>% 
    pivot_wider(names_from = statistic, values_from = value) %>% 
    group_by(iteration, method) %>% 
    summarise(mean_mu = mean(mu),
              mean_sigma2 = mean(sigma2)) %>% 
    ungroup() %>% 
    mutate(sim = type.iter)
  
  sim.plot.df <- rbind(sim.plot.df, combined.df)
}

sim.plot <- ggplot(data = sim.plot.df,
                   mapping = aes(x = mean_mu, y = mean_sigma2, shape = method, colour = method)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(name = "Method",
                     labels = c("ep" = "EP", "laplace" = "Laplace"),
                     values = c("ep" = 19, "laplace" = 17)) +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~sim, scales = "free",
             labeller = as_labeller(c("1" = "n = 200, p.1 = 40, p.2 = 10",
                                      "2" = "n = 200, p.1 = 20, p.2 = 20",
                                      "3" = "n = 200, p.1 = 10, p.2 = 40"))) +
  labs(x = "Mean difference in mu",
       y = "Mean difference in sigma squared",
       shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plots.png", sim.plot, width = 8, height = 5)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-sim-plots.png")

## Benchmarks

bench.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Heteroscedastic/Heteroscedastic-data/Bench-", type.iter, ".RData"))
  p.1 <- ncol(X.1)
  
  load(paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-Laplace.RData"))
  laplace.df <- results.df %>% mutate(method = "laplace")
  load(paste0("Heteroscedastic/Heteroscedastic-results/Bench-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, laplace.df, ep.df) %>% 
    filter(j <= p.1) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_laplace = (mu_laplace - mu_mcmc)/sigma_2_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/sigma_2_mcmc,
           sigma2_laplace = (sigma_2_laplace - sigma_2_mcmc)/sigma_2_mcmc,
           sigma2_ep = (sigma_2_ep - sigma_2_mcmc)/sigma_2_mcmc) %>% 
    select(j, mu_laplace, mu_ep, sigma2_laplace, sigma2_ep) %>% 
    pivot_longer(cols = c(mu_laplace, mu_ep, sigma2_laplace, sigma2_ep)) %>% 
    separate(col = name, into = c("statistic", "method")) %>% 
    pivot_wider(names_from = statistic, values_from = value) %>% 
    mutate(bench = type.iter)
  
  bench.plot.df <- rbind(bench.plot.df, combined.df)
}

bench.plot <- ggplot(data = bench.plot.df,
                     mapping = aes(x = mu, y = sigma2, shape = method, colour = method)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_shape_manual(name = "Method",
                     labels = c("ep" = "EP", "laplace" = "Laplace"),
                     values = c("ep" = 19, "laplace" = 17)) +
  scale_colour_manual(name = "Method",
                      labels = c("ep" = "EP", "laplace" = "Laplace"),
                      values = c("ep" = "black", "laplace" = "darkgray")) +
  facet_wrap(~bench, scales = "free",
             labeller = as_labeller(c("1" = "Food (n = 40, p.1 = 2, p.2 = 2)",
                                      "2" = "Salary (n = 725, p.1 = 6, p.2 = 2)",
                                      "3" = "Sniffer (n = 125, p.1 = 5, p.2 = 5)"))) +
  labs(x = "Difference in mu",
       y = "Difference in sigma squared",
       shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-bench-plots.png", bench.plot, width = 8, height = 5)
plot_crop("Heteroscedastic/Heteroscedastic-plots/Heteroscedastic-bench-plots.png")
