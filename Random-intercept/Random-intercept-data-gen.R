
# Data generation for random intercept linear regression

source("EP-general-auxiliaries.R")
source("Random-intercept/Random-intercept-auxiliaries.R")

set.seed(1)

## Simulation 1

n <- 200
p <- 40
q <- 10
beta <- rep(c(2, -2)/p, p/2)
Z <- as.matrix(bdiag(rep(list(rep(1, n/q)), q)))

for (iteration in 1:num.sim) {
  X <- cbind(1, scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  u <- rnorm(q, sd = exp(kappa.u))
  y <- rnorm(n, X%*%beta + Z%*%u, exp(kappa.y))
  save(X, Z, y, file = paste0("Random-intercept/Random-intercept-data/Sim-1-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
}

## Simulation 2

n <- 200
p <- 20
q <- 20
beta <- rep(c(2, -2)/p, p/2)
Z <- as.matrix(bdiag(rep(list(rep(1, n/q)), q)))

for (iteration in 1:num.sim) {
  X <- cbind(1, scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  u <- rnorm(q, sd = exp(kappa.u))
  y <- rnorm(n, X%*%beta + Z%*%u, exp(kappa.y))
  save(X, Z, y, file = paste0("Random-intercept/Random-intercept-data/Sim-2-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
}

## Simulation 3

n <- 200
p <- 10
q <- 40
beta <- rep(c(2, -2)/p, p/2)
Z <- as.matrix(bdiag(rep(list(rep(1, n/q)), q)))

for (iteration in 1:num.sim) {
  X <- cbind(1, scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
  u <- rnorm(q, sd = exp(kappa.u))
  y <- rnorm(n, X%*%beta + Z%*%u, exp(kappa.y))
  save(X, Z, y, file = paste0("Random-intercept/Random-intercept-data/Sim-3-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
}

## Benchmark 1

load("Benchmark-data/ChickWeight.RData")
X <- as.matrix(cbind(1, createDummyFeatures(ChickWeight[, 4], method = "reference"), scale(ChickWeight[, 2])))
Z <- as.matrix(createDummyFeatures(factor(ChickWeight[, 3], levels = 1:48), method = "1-of-n"))
y <- as.vector(scale(ChickWeight[, 1]))
save(X, Z, y, file = "Random-intercept/Random-intercept-data/Bench-1.RData")

## Benchmark 2

load("Benchmark-data/PEFR.RData")
PEFR <- PEFR %>% arrange(item)
X <- as.matrix(cbind(1, createDummyFeatures(PEFR[, 1], method = "reference")))
Z <- as.matrix(createDummyFeatures(factor(PEFR[, 2]), method = "1-of-n"))
y <- as.vector(scale(PEFR[, 4]))
save(X, Z, y, file = "Random-intercept/Random-intercept-data/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/jsp.RData")
jsp <- jsp %>% filter(school %in% as.character(1:10)) %>% mutate(school = factor(school, levels = 1:10))
X <- as.matrix(cbind(1, createDummyFeatures(jsp[, 3], method = "reference"), scale(jsp[, 7:9])))
Z <- as.matrix(createDummyFeatures(jsp[, 1], method = "1-of-n"))
y <- as.vector(scale(jsp[, 5]))
save(X, Z, y, file = "Random-intercept/Random-intercept-data/Bench-3.RData")
