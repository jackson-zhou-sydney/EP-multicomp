
# General auxiliary functions and variables for expectation propagation

library(tidyverse) # Data manipulation and plotting
library(glmnet)    # Fast ridge regression estimates
library(rstan)     # Robust MCMC and diagnostics
library(EnvStats)  # Compute empirical densities

force.sym <- function(m) {
  # Force matrix to be symmetric
  m.sym <- m
  m.sym[lower.tri(m.sym)] <- t(m.sym)[lower.tri(m.sym)]
  return(m.sym)
}

GI.0 <- function(a, b, c) {
  # Gaussian integral (0th raw moment)
  sqrt(2*pi/a)*exp(-0.5*(c - b^2/(4*a)))
}

GI.1 <- function(a, b, c) {
  # Gaussian integral (1st raw moment)
  sqrt(2*pi/a)*(-b/(2*a))*exp(-0.5*(c - b^2/(4*a)))
}

GI.2 <- function(a, b, c) {
  # Gaussian integral (2nd raw moment)
  sqrt(2*pi/a)*(1/a + b^2/(4*a^2))*exp(-0.5*(c - b^2/(4*a)))
}

TGI.lower.0 <- function(a, b, c) {
  # Truncated Gaussian integral (lower, 0th raw moment)
  exp(log(sqrt(2*pi/a)) - 0.5*(c - b^2/(4*a)) + pnorm(0, -b/(2*a), sqrt(1/a), log.p = TRUE))
}

TGI.lower.1 <- function(a, b, c) {
  # Truncated Gaussian integral (lower, 1st raw moment)
  sqrt(2*pi/a)*(-(b/(2*a))*exp(-0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)), log.p = TRUE)) - (1/sqrt(a))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)), log = TRUE)))
}

TGI.lower.2 <- function(a, b, c) {
  # Truncated Gaussian integral (lower, 2nd raw moment)
  sqrt(2*pi/a)*((b^2/(4*a) + 1)*exp(-0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)), log.p = TRUE)) + (b/(2*sqrt(a)))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)), log = TRUE)))/a
}

TGI.upper.0 <- function(a, b, c) {
  # Truncated Gaussian integral (upper, 0th raw moment)
  exp(log(sqrt(2*pi/a)) - 0.5*(c - b^2/(4*a)) + pnorm(0, -b/(2*a), sqrt(1/a), lower.tail = FALSE, log.p = TRUE))
}

TGI.upper.1 <- function(a, b, c) {
  # Truncated Gaussian integral (upper, 1st raw moment)
  sqrt(2*pi/a)*(-(b/(2*a))*exp(-0.5*(c - b^2/(4*a)) + pnorm(-b/(2*sqrt(a)), log.p = TRUE)) + (1/sqrt(a))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)), log = TRUE)))
}

TGI.upper.2 <- function(a, b, c) {
  # Truncated Gaussian integral (upper, 2nd raw moment)
  sqrt(2*pi/a)*((b^2/(4*a) + 1)*exp(-0.5*(c - b^2/(4*a)) + pnorm(-b/(2*sqrt(a)), log.p = TRUE)) - (b/(2*sqrt(a)))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)), log = TRUE)))/a
}
