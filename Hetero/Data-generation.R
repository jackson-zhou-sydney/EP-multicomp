
# Data generation for heteroscedastic linear regression

source("General-auxiliaries.R")
source("Hetero/Auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.sim) {
  n <- sim.settings[[type.iter]][["n"]]
  p.1 <- sim.settings[[type.iter]][["p.1"]]
  p.2 <- sim.settings[[type.iter]][["p.2"]]
  beta.1 <- rep(c(2, -2)/p.1, p.1/2)
  beta.2 <- rep(c(2, -2)/p.2, p.2/2)
  
  for (iteration in 1:num.sim.iter) {
    X.1 <- cbind(1, scale(matrix(data = rnorm(n = n*(p.1 - 1)), nrow = n)))
    X.2 <- cbind(1, scale(matrix(data = rnorm(n = n*(p.2 - 1)), nrow = n)))
    y <- as.vector(scale(rnorm(n, X.1%*%beta.1, exp(X.2%*%beta.2))))
    save(X.1, X.2, y, file = paste0("Hetero/Data/Simulations/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
  }
}

## Benchmark 1

foodexp <- as.matrix(read_dta("Benchmark-data/foodexp.dta"))
X.1 <- cbind(1, scale(foodexp[, 2]))
X.2 <- cbind(1, scale(foodexp[, 2]))
y <- as.vector(scale(foodexp[, 1]))
save(X.1, X.2, y, file = "Hetero/Data/Benchmarks/Bench-1.RData")

## Benchmark 2

salary <- as.matrix(read_dta("Benchmark-data/salary.dta"))
X.1 <- cbind(1, salary[, 13], unname(scale(salary[, 3:6])))
X.2 <- cbind(1, salary[, 13])
y <- as.vector(scale(salary[, 1]))
save(X.1, X.2, y, file = "Hetero/Data/Benchmarks/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/sniffer.RData")
X.1 <- cbind(1, scale(unname(sniffer[, -5])))
X.2 <- cbind(1, scale(unname(sniffer[, -5])))
y <- as.vector(scale(sniffer[, 5]))
save(X.1, X.2, y, file = "Hetero/Data/Benchmarks/Bench-3.RData")

## Big data

energy <- read.csv("Benchmark-data/energydata_complete.csv")
energy_sub <- energy[, -c(1, 3, 28, 29)]
energy_squared <- energy_sub %>% select(-Appliances) %>% mutate_all(function(x) x^2) 
colnames(energy_squared) <- paste0(colnames(energy_squared), "_sqr")
energy_all <- cbind(energy_sub, energy_squared)

X.1.all <- unname(model.matrix(Appliances ~ .^2, data = energy_all))
attr(X.1.all, "assign") <- NULL
X.2.all <- unname(as.matrix(cbind(1, energy_sub[, -1])))
y.all <- as.vector(energy_all[, 1])

train.id <- sample(1:nrow(energy), 3000)

for (j in 2:ncol(X.1.all)) {
  X.1.all[, j] <- X.1.all[, j] - mean(X.1.all[train.id, j])
  X.1.all[, j] <- X.1.all[, j]/sd(X.1.all[train.id, j])
}

for (j in 2:ncol(X.2.all)) {
  X.2.all[, j] <- X.2.all[, j] - mean(X.2.all[train.id, j])
  X.2.all[, j] <- X.2.all[, j]/sd(X.2.all[train.id, j])
}

y.all <- y.all - mean(y.all[train.id])
y.all <- y.all/sd(y.all[train.id])

X.1 <- X.1.all[train.id, ]
X.2 <- X.2.all[train.id, ]
y <- y.all[train.id]

save(X.1, X.2, y, file = "Hetero/Data/Big/Big.RData")

for (i in 1:8) {
  test.id <- sample(setdiff(1:nrow(energy), train.id))[1:ceiling(train.size*3000)]
  
  X.1.test <- X.1.all[test.id, ]
  X.2.test <- X.2.all[test.id, ]
  y.test <- y.all[test.id]
  
  save(X.1.test, X.2.test, y.test, file = "Hetero/Data/Big/Big-test-", str_pad(i, 2, pad = "0"), ".RData")
}
