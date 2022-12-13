
# General auxiliary functions and variables for expectation propagation

library(tidyverse)   # Data manipulation and plotting
library(rstan)       # Robust MCMC and diagnostics
library(EnvStats)    # Compute empirical densities
library(glmnet)      # Fast ridge regression estimates
library(quantreg)    # Quantile regression estimates
library(lme4)        # Fast linear mixed model estimates
library(stringr)     # String manipulation
library(knitr)       # Crop plots
library(haven)       # Load in DTA files
library(mlr)         # Generate dummy columns
library(mvtnorm)     # Sample from multivariate Gaussians
library(nbpMatching) # Non-bipartite matching
library(pracma)      # Trapezoidal integration

num.each.type <- 3
num.sim <- 50
mcmc.chains <- 1
mcmc.iter <- 50000
mcmc.warmup <- 5000
mcmc.short.iter <- 2000
mcmc.short.warmup <- 200
match.size <- 1000
match.reps <- 10
sd.multiple <- 5
total.grid.points <- 1000

err <- function(e) {
  # Return NA on error
  return(NA)
}

sym <- function(m) {
  # Force matrix to be symmetric
  m.sym <- m
  m.sym[lower.tri(m.sym)] <- t(m.sym)[lower.tri(m.sym)]
  return(m.sym)
}

nbp.match.pairs <- function(X, Y) {
  # Non-bipartite match pairs
  X.size <- nrow(X)
  Y.size <- nrow(Y)
  df <- data.frame(rbind(X, Y))
  df.rank <- apply(df, 2, function(s) rank(s, ties.method = "first"))
  df.rank <- data.frame(cbind(1:nrow(df.rank), df.rank))
  nbp.res <- nonbimatch(distancematrix(gendistance(df.rank, 1)))
  return(sum((c(rep(1, X.size), rep(0, Y.size))[nbp.res$matches[1:X.size, 4]]) != 1))
}

block.mean <- function(m, l, na.rm) {
  # Block-wise mean of matrix
  ret.mat <- matrix(nrow = length(l), ncol = length(l))
  for (i in 1:length(l)) {
    for (j in 1:length(l)) {
      ret.mat[i, j] <- mean(m[l[[i]], l[[j]]], na.rm = na.rm)
    }
  }
  return(ret.mat)
}

GI.0 <- function(x) {
  # Gaussian integral (0th raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*exp(-0.5*(c - b^2/(4*a)))
}

GI.1 <- function(x) {
  # Gaussian integral (1st raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*(-b/(2*a))*exp(-0.5*(c - b^2/(4*a)))
}

GI.2 <- function(x) {
  # Gaussian integral (2nd raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*(1/a + b^2/(4*a^2))*exp(-0.5*(c - b^2/(4*a)))
}

TGI.minus.0 <- function(x, y) {
  # Truncated Gaussian integral (lower, 0th raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  exp(log(sqrt(2*pi/a)) - 0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, log.p = T))
}

TGI.minus.1 <- function(x, y) {
  # Truncated Gaussian integral (lower, 1st raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*(-(b/(2*a))*exp(-0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, log.p = T)) - (1/sqrt(a))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))
}

TGI.minus.2 <- function(x, y) {
  # Truncated Gaussian integral (lower, 2nd raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*((b^2/(4*a) + 1)*exp(-0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, log.p = T)) + (b/(2*sqrt(a)) - sqrt(a)*y)*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))/a
}

TGI.plus.0 <- function(x, y) {
  # Truncated Gaussian integral (upper, 0th raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  exp(log(sqrt(2*pi/a)) - 0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, lower.tail = F, log.p = T))
}

TGI.plus.1 <- function(x, y) {
  # Truncated Gaussian integral (upper, 1st raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*(-(b/(2*a))*exp(-0.5*(c - b^2/(4*a)) + pnorm(-b/(2*sqrt(a)) - sqrt(a)*y, log.p = T)) + (1/sqrt(a))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))
}

TGI.plus.2 <- function(x, y) {
  # Truncated Gaussian integral (upper, 2nd raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*((b^2/(4*a) + 1)*exp(-0.5*(c - b^2/(4*a)) + pnorm(-b/(2*sqrt(a)) - sqrt(a)*y, log.p = T)) - (b/(2*sqrt(a)) - sqrt(a)*y)*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))/a
}
