
# Data generation for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  n <- sim.settings[[type.iter]][["n"]]
  p <- sim.settings[[type.iter]][["p"]]
  beta <- rep(c(2, -2)/p, p/2)
  beta[1:p/2] <- 0
  
  for (iteration in 1:num.sim) {
    X <- cbind(1, scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
    y <- rnorm(n, X%*%beta, exp(kappa))
    save(X, y, file = paste0("Lasso/Lasso-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
  }
}

## Benchmark 1

load("Benchmark-data/efron2004.Rdata")
X <- cbind(1, scale(efron2004$x))
y <- as.vector(efron2004$y)
save(X, y, file = "Lasso/Lasso-data/Bench-1.RData")

## Benchmark 2

load("Benchmark-data/Prostate.RData")
X <- cbind(1, scale(Prostate[, -9]))
y <- Prostate[, 9]
save(X, y, file = "Lasso/Lasso-data/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/eyedata.RData")
X <- cbind(1, scale(unname(x)))
y <- y
save(X, y, file = "Lasso/Lasso-data/Bench-3.RData")
