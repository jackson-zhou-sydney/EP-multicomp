
# Pathfinder debugging

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

seed <- 1
set.seed(seed)

library(cmdstanr)
mcmc <- cmdstan_model("Hetero/Methods/MCMC.stan", stanc_options = list("O1"))

type.iter <- 1
iteration <- 1

load(paste0("Hetero/Data/Simulations/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
n <- nrow(X.1)
p.1 <- ncol(X.1)
p.2 <- ncol(X.2)

mu.theta <- rep(0, p.1 + p.2)
Sigma.theta <- diag(c(rep(sigma.2.beta.1, p.1), rep(sigma.2.beta.2, p.2)))

train.ind <- sample(1:n)[1:ceiling(train.size*n)]
X.1.train <- X.1[train.ind, , drop = F]
X.2.train <- X.2[train.ind, , drop = F]
y.train <- y[train.ind]
n.train <- nrow(X.1.train)

X.1.test <- X.1[-train.ind, , drop = F]
X.2.test <- X.2[-train.ind, , drop = F]
y.test <- y[-train.ind]
n.test <- nrow(X.1.test)

mcmc.rstan <- rstan::stan_model("Hetero/Methods/MCMC.stan")

load("Hetero/Results/Simulations-conv-table.RData")
mcmc.test.iter <- sim.r.hat.table %>% pull(mcmc_iter) %>% unique() %>% sort()
iters <- sim.r.hat.table %>% filter(sim == type.iter) %>% filter(mean_max_r_hat > r.hat.tol) %>% pull(mcmc_iter)
ind <- if(length(iters) != 0) which(mcmc.test.iter == max(iters)) + 1 else 1
mcmc.s.iter <- mcmc.test.iter[min(ind, length(mcmc.test.iter))]

start.time <- proc.time()
Rprof(tf <- "rprof.log", memory.profiling = T, interval = 0.005)

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
                                N_sam = round(mcmc.s.iter/num.cores),
                                init_bound = 0.1)

gvb.samples <- t(Imp_Resam_WR(opath, n_sam = mcmc.s.iter, seed = seed))
gvb.mu <- colMeans(gvb.samples)
gvb.Sigma <- var(gvb.samples)

total.time <- proc.time() - start.time
Rprof(NULL)

print(total.time)
summaryRprof(tf, memory = "stats")
