
# Data generation for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

set.seed(1)

## Simulation 1

n <- 200
p <- 40
beta <- rep(c(2, -2)/p, p/2)
beta[1:p/2] <- 0

for (iteration in 1:num.sim) {
  X <- cbind(rep(1, n), scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  y <- rnorm(n, X%*%beta, sigma)
  save(X, y, file = paste0("Lasso/Lasso-data/Sim-1-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
}

## Simulation 2

n <- 40
p <- 40
beta <- rep(c(2, -2)/p, p/2)
beta[1:p/2] <- 0

for (iteration in 1:num.sim) {
  X <- cbind(rep(1, n), scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  y <- rnorm(n, X%*%beta, sigma)
  save(X, y, file = paste0("Lasso/Lasso-data/Sim-2-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
}

## Simulation 3

n <- 10
p <- 40
beta <- rep(c(2, -2)/p, p/2)
beta[1:p/2] <- 0

for (iteration in 1:num.sim) {
  X <- cbind(rep(1, n), scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  y <- rnorm(n, X%*%beta, sigma)
  save(X, y, file = paste0("Lasso/Lasso-data/Sim-3-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
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
