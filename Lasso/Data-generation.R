
# Data generation for lasso linear regression

source("General-auxiliaries.R")
source("Lasso/Auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.sim) {
  n <- sim.settings[[type.iter]][["n"]]
  p <- sim.settings[[type.iter]][["p"]]
  beta <- rep(c(2, -2)/p, p/2)
  beta[1:p/2] <- 0
  
  for (iteration in 1:num.sim.iter) {
    X <- cbind(1, scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
    y <- rnorm(n, X%*%beta, exp(kappa))
    save(X, y, file = paste0("Lasso/Data/Simulations/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
  }
}

## Benchmark 1

load("Benchmark-data/efron2004.Rdata")
X <- cbind(1, scale(efron2004$x))
y <- as.vector(efron2004$y)
save(X, y, file = "Lasso/Data/Benchmarks/Bench-1.RData")

## Benchmark 2

load("Benchmark-data/Prostate.RData")
X <- cbind(1, scale(Prostate[, -9]))
y <- Prostate[, 9]
save(X, y, file = "Lasso/Data/Benchmarks/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/eyedata.RData")
X <- cbind(1, scale(unname(x)))
y <- y
save(X, y, file = "Lasso/Data/Benchmarks/Bench-3.RData")

## Benchmark 4

energy <- read.csv("Benchmark-data/energydata_complete.csv")
energy_cleaned <- as.data.frame(scale(energy[, -c(1, 3, 28, 29)]))

energy_squared <- energy_cleaned %>% 
  mutate_all(function(x) x^2) %>% 
  select(-Appliances)
colnames(energy_squared) <- paste0(colnames(energy_squared), "_sqr")

energy_all <- cbind(energy_cleaned, energy_squared)

X <- unname(model.matrix(Appliances ~ .^2, data = energy_all))
attr(X, "assign") <- NULL
X[, 2:1177] <- scale(X[, 2:1177])
y <- energy_all[, 1]
save(X, y, file = "Lasso/Data/Benchmarks/Bench-4.RData")
