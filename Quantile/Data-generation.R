
# Data generation for quantile regression

source("General-auxiliaries.R")
source("Quantile/Auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.sim) {
  n <- as.numeric(sim.settings[[type.iter]][["n"]])
  p <- as.numeric(sim.settings[[type.iter]][["p"]])
  dist <- sim.settings[[type.iter]][["dist"]]
  beta <- rep(c(2, -2)/p, p/2)
  
  for (iteration in 1:num.sim.iter) {
    X <- cbind(1, scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
    
    if (dist == "normal") {
      y <- as.vector(scale(rnorm(n, X%*%beta, sigma)))
    } else if (dist == "poisson") {
      y <- as.vector(scale(rpois(n, exp(X%*%beta))))
    } else if (dist == "binomial") {
      y <- as.vector(scale(rbinom(n, n.binom, pnorm(X%*%beta))))
    }
    
    save(X, y, file = paste0("Quantile/Data/Simulations/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
  }
}

## Benchmark 1

load("Benchmark-data/ImmunogG.RData")
X <- cbind(1, scale(ImmunogG$Age))
y <- as.vector(scale(ImmunogG$IgG))
save(X, y, file = "Quantile/Data/Benchmarks/Bench-1.RData")

## Benchmark 2

load("Benchmark-data/engel.RData")
X <- cbind(1, scale(engel[, 1]))
y <- as.vector(scale(engel[, 2]))
save(X, y, file = "Quantile/Data/Benchmarks/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/stackloss.RData")
X <- cbind(1, scale(stackloss[, 1:3]))
y <- as.vector(scale(stackloss[, 4]))
save(X, y, file = "Quantile/Data/Benchmarks/Bench-3.RData")

## Benchmark 4

energy <- read.csv("Benchmark-data/energydata_complete.csv")
energy_filtered <- energy[sample(1:nrow(energy), 6000), ]
energy_cleaned <- as.data.frame(scale(energy_filtered[, -c(1, 3, 28, 29)]))

energy_squared <- energy_cleaned %>% 
  mutate_all(function(x) x^2) %>% 
  select(-Appliances)
colnames(energy_squared) <- paste0(colnames(energy_squared), "_sqr")

energy_all <- cbind(energy_cleaned, energy_squared)

X <- unname(model.matrix(Appliances ~ .^2, data = energy_all))
attr(X, "assign") <- NULL
X[, 2:1177] <- scale(X[, 2:1177])
y <- energy_all[, 1]
save(X, y, file = "Quantile/Data/Benchmarks/Bench-4.RData")
