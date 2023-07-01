
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
    y <- as.vector(scale(rnorm(n, X%*%beta, exp(kappa))))
    save(X, y, file = paste0("Lasso/Data/Simulations/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
  }
}

## Benchmark 1

load("Benchmark-data/efron2004.Rdata")
X <- cbind(1, scale(efron2004$x))
y <- as.vector(scale(efron2004$y))
save(X, y, file = "Lasso/Data/Benchmarks/Bench-1.RData")

## Benchmark 2

load("Benchmark-data/Prostate.RData")
X <- cbind(1, scale(Prostate[, -9]))
y <- as.vector(scale(Prostate[, 9]))
save(X, y, file = "Lasso/Data/Benchmarks/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/eyedata.RData")
X <- cbind(1, scale(unname(x)))
y <- as.vector(scale(y))
save(X, y, file = "Lasso/Data/Benchmarks/Bench-3.RData")

## Big data

energy <- read.csv("Benchmark-data/energydata_complete.csv")
energy_sub <- energy[, -c(1, 3, 28, 29)]
energy_squared <- energy_sub %>% select(-Appliances) %>% mutate_all(function(x) x^2) 
colnames(energy_squared) <- paste0(colnames(energy_squared), "_sqr")
energy_all <- cbind(energy_sub, energy_squared)

X.all <- unname(model.matrix(Appliances ~ .^2, data = energy_all))
attr(X.all, "assign") <- NULL
y.all <- as.vector(energy_all[, 1])

train.id <- sample(1:nrow(energy), round(train.size*nrow(energy)))

for (j in 2:ncol(X.all)) {
  X.all[, j] <- X.all[, j] - mean(X.all[train.id, j])
  X.all[, j] <- X.all[, j]/sd(X.all[train.id, j])
}

y.all <- y.all - mean(y.all[train.id])
y.all <- y.all/sd(y.all[train.id])

X <- X.all[train.id, ]
y <- y.all[train.id]

X.test <- X.all[-train.id, ]
y.test <- y.all[-train.id]

save(X, y, file = "Lasso/Data/Big/Big.RData")
save(X.test, y.test, file = "Lasso/Data/Big/Big-test.RData")
