
# Data generation for heteroscedastic linear regression

source("EP-general-auxiliaries.R")
source("Heteroscedastic/Heteroscedastic-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  n <- sim.settings[[type.iter]][["n"]]
  p.1 <- sim.settings[[type.iter]][["p.1"]]
  p.2 <- sim.settings[[type.iter]][["p.2"]]
  beta.1 <- rep(c(2, -2)/p.1, p.1/2)
  beta.2 <- rep(c(2, -2)/p.2, p.2/2)
  
  for (iteration in 1:num.sim) {
    X.1 <- cbind(1, scale(matrix(data = rnorm(n = n*(p.1 - 1)), nrow = n)))
    X.2 <- cbind(1, scale(matrix(data = rnorm(n = n*(p.2 - 1)), nrow = n)))
    y <- as.vector(scale(rnorm(n, X.1%*%beta.1, sqrt(exp(X.2%*%beta.2)))))
    save(X.1, X.2, y, file = paste0("Heteroscedastic/Heteroscedastic-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
  }
}

## Benchmark 1

foodexp <- as.matrix(read_dta("Benchmark-data/foodexp.dta"))
X.1 <- cbind(1, scale(foodexp[, 2]))
X.2 <- cbind(1, scale(foodexp[, 2]))
y <- as.vector(scale(foodexp[, 1]))
save(X.1, X.2, y, file = "Heteroscedastic/Heteroscedastic-data/Bench-1.RData")

## Benchmark 2

salary <- as.matrix(read_dta("Benchmark-data/salary.dta"))
X.1 <- cbind(1, salary[, 13], unname(scale(salary[, 3:6])))
X.2 <- cbind(1, salary[, 13])
y <- as.vector(scale(salary[, 1]))
save(X.1, X.2, y, file = "Heteroscedastic/Heteroscedastic-data/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/sniffer.RData")
X.1 <- cbind(1, scale(unname(sniffer[, -5])))
X.2 <- cbind(1, scale(unname(sniffer[, -5])))
y <- as.vector(scale(sniffer[, 5]))
save(X.1, X.2, y, file = "Heteroscedastic/Heteroscedastic-data/Bench-3.RData")
