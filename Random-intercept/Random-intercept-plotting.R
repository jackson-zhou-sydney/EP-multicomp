
# Plots for random intercept linear regression

source("EP-general-auxiliaries.R")
source("Random-intercept/Random-intercept-auxiliaries.R")

## Simulations

sim.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Random-intercept/Random-intercept-data/Sim-", type.iter, "-iter-01.RData"))
  p <- ncol(X)
  
  load(paste0("Random-intercept/Random-intercept-results/Sim-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Random-intercept/Random-intercept-results/Sim-", type.iter, "-res-Laplace.RData"))
  laplace.df <- results.df %>% mutate(method = "laplace")
  load(paste0("Random-intercept/Random-intercept-results/Sim-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, laplace.df, ep.df) %>% 
    filter(j <= p) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_laplace = (mu_laplace - mu_mcmc)/mu_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/mu_mcmc,
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
             labeller = as_labeller(c("1" = "n = 200, p = 40, q = 10",
                                      "2" = "n = 200, p = 20, q = 20",
                                      "3" = "n = 200, p = 10, q = 40"))) +
  labs(x = "Mean difference in mu",
       y = "Mean difference in sigma squared",
       shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Random-intercept/Random-intercept-plots/Random-intercept-sim-plots.png", sim.plot, width = 8, height = 5)
plot_crop("Random-intercept/Random-intercept-plots/Random-intercept-sim-plots.png")

## Benchmarks

bench.plot.df <- data.frame()

for (type.iter in 1:num.each.type) {
  load(paste0("Random-intercept/Random-intercept-data/Bench-", type.iter, ".RData"))
  p <- ncol(X)
  
  load(paste0("Random-intercept/Random-intercept-results/Bench-", type.iter, "-res-MCMC.RData"))
  mcmc.df <- results.df %>% mutate(method = "mcmc")
  load(paste0("Random-intercept/Random-intercept-results/Bench-", type.iter, "-res-Laplace.RData"))
  laplace.df <- results.df %>% mutate(method = "laplace")
  load(paste0("Random-intercept/Random-intercept-results/Bench-", type.iter, "-res-EP.RData"))
  ep.df <- results.df %>% mutate(method = "ep")
  
  combined.df <- rbind(mcmc.df, laplace.df, ep.df) %>% 
    filter(j <= p) %>% 
    pivot_wider(names_from = method, values_from = c(mu, sigma_2)) %>% 
    mutate(mu_laplace = (mu_laplace - mu_mcmc)/mu_mcmc,
           mu_ep = (mu_ep - mu_mcmc)/mu_mcmc,
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
             labeller = as_labeller(c("1" = "Chick (n = 578, p = 5, q = 50)",
                                      "2" = "PEFR (n = 68, p = 2, q = 17)",
                                      "3" = "JSP (n = 606, p = 5, q = 10)"))) +
  labs(x = "Difference in mu",
       y = "Difference in sigma squared",
       shape = "Method") +
  theme_bw() +
  theme(legend.position = "top")

ggsave("Random-intercept/Random-intercept-plots/Random-intercept-bench-plots.png", bench.plot, width = 8, height = 5)
plot_crop("Random-intercept/Random-intercept-plots/Random-intercept-bench-plots.png")
