
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

## Options
options(rcpp.cache.dir = "Rcpp-cache/")

## General settings
num.cores <- 10
num.sim.iter <- 5
num.sim <- 3
num.bench <- 4

## Marginal L1 evaluation
sd.multiple <- 5
total.grid.points <- 1024

## MMD and LPPD evaluation
eval.size <- 500
train.size <- 0.9

## R hat evaluation
r.hat.tol <- 1.05
warmup.mult <- 0.1
long.mult <- 4
mcmc.check.iter <- c(50, 100, 200, 500, 1000, 2000)
mcmc.check.warmup <- warmup.mult*mcmc.check.iter
