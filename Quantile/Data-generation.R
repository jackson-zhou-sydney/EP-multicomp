
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

## Big data

energy <- read.csv("Benchmark-data/energydata_complete.csv")
energy_sub <- energy[, -c(1, 3, 28, 29)]

X.all <- unname(model.matrix(Appliances ~ .^2, data = energy_sub))
attr(X.all, "assign") <- NULL
y.all <- as.vector(energy_sub[, 1])

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

save(X, y, file = "Quantile/Data/Big/Big.RData")
save(X.test, y.test, file = "Quantile/Data/Big/Big-test.RData")
