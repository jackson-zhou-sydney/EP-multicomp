
# Pathfinder debugging

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

seed <- 1
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Hetero/Methods/MCMC.stan", stanc_options = list("O1"))

load("Hetero/Data/Big/Big.RData")
n <- nrow(X.1)
p.1 <- ncol(X.1)
p.2 <- ncol(X.2)

mu.theta <- rep(0, p.1 + p.2)
Sigma.theta <- diag(c(rep(sigma.2.beta.1, p.1), rep(sigma.2.beta.2, p.2)))

load(paste0("Hetero/Results/Big-MCMC-results-", str_pad(seed + 2, 2, pad = "0"), ".Rdata"))
ind <- which(big.test.iter == r.hat.df %>% filter(max_r_hat > r.hat.tol) %>% pull(mcmc_iter) %>% max()) + 1
mcmc.s.iter <- big.test.iter[min(ind, length(big.test.iter))]

opath <- opt_path_parallel(seed_init = (seed - 1)*length(num.cores) + 1:num.cores, 
                           seed_list = (seed - 1)*length(num.cores) + 1:num.cores, 
                           mc.cores = num.cores, 
                           fn = function(theta) nljl(theta, X.1, X.2, y, Sigma.theta, mu.theta),
                           gr = function(theta) nljl_grad(theta, X.1, X.2, y, Sigma.theta, mu.theta),
                           D = p.1 + p.2,
                           N_sam = round(mcmc.s.iter/num.cores),
                           init_bound = 0.1)

gvb.samples <- t(Imp_Resam_WR(opath, n_sam = mcmc.s.iter, seed = seed))
gvb.mu <- colMeans(gvb.samples)
gvb.Sigma <- var(gvb.samples)

print("Done")
