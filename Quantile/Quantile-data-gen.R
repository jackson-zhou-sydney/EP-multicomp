
# Data generation for quantile regression

source("EP-general-auxiliaries.R")
source("Quantile/Quantile-auxiliaries.R")

set.seed(1)

## Simulations

for (type.iter in 1:num.each.type) {
  n <- as.numeric(sim.settings[[type.iter]][["n"]])
  p <- as.numeric(sim.settings[[type.iter]][["p"]])
  dist <- sim.settings[[type.iter]][["dist"]]
  beta <- rep(c(2, -2)/p, p/2)
  
  for (iteration in 1:num.sim) {
    X <- cbind(1, scale(matrix(data = rnorm(n = n*(p - 1)), nrow = n)))
    
    if (dist == "normal") {
      y <- rnorm(n, X%*%beta, sigma)
    } else if (dist == "poisson") {
      y <- rpois(n, exp(X%*%beta))
    } else if (dist == "binomial") {
      y <- rbinom(n, n.binom, pnorm(X%*%beta))
    }
    
    save(X, y, file = paste0("Quantile/Quantile-data/Sim-", type.iter, "-iter-", str_pad(iteration, 2, pad = "0"), ".RData"))
  }
}

## Benchmark 1

load("Benchmark-data/ImmunogG.RData")
X <- cbind(1, scale(ImmunogG$Age))
y <- ImmunogG$IgG
save(X, y, file = "Quantile/Quantile-data/Bench-1.RData")

## Benchmark 2

load("Benchmark-data/engel.RData")
X <- cbind(1, scale(engel[, 1]))
y <- as.vector(scale(engel[, 2]))
save(X, y, file = "Quantile/Quantile-data/Bench-2.RData")

## Benchmark 3

load("Benchmark-data/stackloss.RData")
X <- cbind(1, scale(stackloss[, 1:3]))
y <- as.vector(scale(stackloss[, 4]))
save(X, y, file = "Quantile/Quantile-data/Bench-3.RData")

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
y <- energy_all$Appliances
save(X, y, file = "Quantile/Quantile-data/Bench-4.RData")
