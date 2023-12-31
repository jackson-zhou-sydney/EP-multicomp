
# General auxiliary functions and variables for expectation propagation

library(methods)       # Standard methods in R
library(tidyverse)     # Data manipulation and plotting
library(EnvStats)      # Compute empirical densities
library(glmnet)        # Fast ridge regression estimates
library(quantreg)      # Quantile regression estimates
library(lme4)          # Fast linear mixed model estimates
library(stringr)       # String manipulation
library(knitr)         # Crop plots
library(haven)         # Load in DTA files
library(mlr)           # Generate dummy columns
library(mvtnorm)       # Sample from multivariate Gaussians
library(nbpMatching)   # Non-bipartite matching
library(pracma)        # Trapezoidal integration
library(kernlab)       # Maximum mean discrepancy
library(Rcpp)          # C++ in R
library(RcppArmadillo) # Efficient linear algebra for Rcpp
library(RcppEigen)     # Alternate linear algebra for Rcpp
library(RcppNumerical) # Optimisation in Rcpp
library(calculus)      # Numerical derivatives

## Options
options(rcpp.cache.dir = "Rcpp-cache/")

## General settings
num.cores <- 10
num.sim.iter <- 5
num.sim <- 3
num.bench <- 3
log.lb <- 0.00001
table.dp <- 2
l1.dp <- 1
train.size <- 0.875

## Marginal L1 evaluation
sd.multiple <- 5
total.grid.points <- 1024

## MMD and LPPD evaluation
eval.size <- 500

## R hat evaluation
r.hat.tol <- 1.1
warmup.mult <- 0.1
mcmc.g.iter <- 10000
mcmc.g.warmup <- warmup.mult*mcmc.g.iter
big.test.iter <- c(100, 200, 400, 1000, 2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000, 18000, 20000)
big.mcmc.g.iter <- 20000
big.mcmc.g.warmup <- warmup.mult*big.mcmc.g.iter
max.tree.depth <- 10
