
# Data generation for lasso linear regression

source("EP-general-auxiliaries.R")
source("Lasso/Lasso-auxiliaries.R")

set.seed(1)
sigma <- 0.1

## Simulation 1

n <- 200
p <- 40
beta <- rep(c(2, -2)/p, p/2)
beta[1:p/2] <- 0

for (iteration in 1:num.sim) {
  X <- cbind(rep(1, n), scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  y <- rnorm(n, X%*%beta, sigma)
  write.table(cbind(X, y), file = paste0("Lasso/Lasso-data/Sim-1-iter-", str_pad(iteration, 2, pad = "0"), ".csv"),
              row.names = F, col.names = F, sep = ",")
}

## Simulation 2

n <- 40
p <- 40
beta <- rep(c(2, -2)/p, p/2)
beta[1:p/2] <- 0

for (iteration in 1:num.sim) {
  X <- cbind(rep(1, n), scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  y <- rnorm(n, X%*%beta, sigma)
  write.table(cbind(X, y), file = paste0("Lasso/Lasso-data/Sim-2-iter-", str_pad(iteration, 2, pad = "0"), ".csv"),
              row.names = F, col.names = F, sep = ",")
}

## Simulation 3

n <- 10
p <- 40
beta <- rep(c(2, -2)/p, p/2)
beta[1:p/2] <- 0

for (iteration in 1:num.sim) {
  X <- cbind(rep(1, n), scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  y <- rnorm(n, X%*%beta, sigma)
  write.table(cbind(X, y), file = paste0("Lasso/Lasso-data/Sim-3-iter-", str_pad(iteration, 2, pad = "0"), ".csv"),
              row.names = F, col.names = F, sep = ",")
}

## Benchmark 1

load("Benchmark-data/efron2004.Rdata")
X <- cbind(1, scale(efron2004$x))
y <- as.vector(efron2004$y)
write.table(cbind(X, y), file = "Lasso/Lasso-data/Bench-1.csv",
            row.names = F, col.names = F, sep = ",")

## Benchmark 2

load("Benchmark-data/Prostate.RData")
X <- cbind(1, scale(Prostate[, -9]))
y <- Prostate[, 9]
write.table(cbind(X, y), file = "Lasso/Lasso-data/Bench-2.csv",
            row.names = F, col.names = F, sep = ",")

## Benchmark 3

load("Benchmark-data/eyedata.RData")
X <- cbind(1, scale(unname(x)))
y <- y
write.table(cbind(X, y), file = "Lasso/Lasso-data/Bench-3.csv",
            row.names = F, col.names = F, sep = ",")
