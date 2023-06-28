
source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

type.iter <- 1
iteration <- 1
seed <- 1
set.seed(1)

library(cmdstanr)
mcmc <- cmdstan_model("Hetero/Methods/MCMC.stan", stanc_options = list("O1"))

load(paste0("Hetero/Data/Simulations/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
n <- nrow(X.1)
p.1 <- ncol(X.1)
p.2 <- ncol(X.2)

mu.theta <- rep(0, p.1 + p.2)
Sigma.theta <- diag(c(rep(sigma.2.beta.1, p.1), rep(sigma.2.beta.2, p.2)))

mcmc.rstan <- rstan::stan_model("Hetero/Methods/MCMC.stan")
mcmc.s.iter <- 2000
# num.cores <- 1

seed <- 1
opath <- opt_path_stan_parallel(seed_init = (seed - 1)*length(num.cores) + 1:num.cores, 
                                seed_list = (seed - 1)*length(num.cores) + 1:num.cores, 
                                mc.cores = num.cores, 
                                model = mcmc.rstan,
                                data = list(N = n,
                                            p_1 = p.1,
                                            p_2 = p.2,
                                            X_1 = X.1,
                                            X_2 = X.2,
                                            y = y,
                                            mu_theta = mu.theta,
                                            Sigma_theta = Sigma.theta),
                                N_sam = round(mcmc.s.iter/num.cores))

# gvb.samples <- t(Imp_Resam_WR(opath, n_sam = mcmc.s.iter, seed = seed))
# gvb.mu <- colMeans(gvb.samples)
# gvb.Sigma <- var(gvb.samples)
