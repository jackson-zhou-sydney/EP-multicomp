
# Auxiliary functions and variables for lasso linear regression

sourceCpp("Lasso/Methods/EP.cpp")
sourceCpp("Lasso/Methods/EP-2D.cpp")
sourceCpp("Lasso/Methods/MFVB.cpp")

mu.kappa <- 0
sigma.2.kappa <- 0.01
lambda <- 0.5
kappa <- -1

sim.settings <- list(c(n = 200, p = 40),
                     c(n = 40, p = 40),
                     c(n = 10, p = 40))

sim.labels <- c("1" = "n = 200, p = 40",
                "2" = "n = 40, p = 40",
                "3" = "n = 10, p = 40")

bench.settings <- list(c(n = 442, p = 11),
                       c(n = 97, p = 9),
                       c(n = 120, p = 201),
                       c(n = 17268, p = 667))

bench.labels <- c("1" = "Diabetes (n = 442, p = 11)",
                  "2" = "Prostate (n = 97, p = 9)",
                  "3" = "Eye (n = 120, p = 201)",
                  "4" = "Energy (n = 17268, p = 667)")

point.likelihood <- function(theta, x, y) {
  # Likelihood evaluated at a point
  p <- length(x)
  
  beta <- theta[1:p]
  kappa <- theta[p + 1]
  
  as.numeric(exp(-kappa - (y - t(x)%*%beta)^2/(2*exp(2*kappa)))/sqrt(2*pi))
}

lppd <- function(X, y, S) {
  # Compute the log pointwise predictive density
  # S is a matrix of posterior samples
  total <- 0
  n <- nrow(X)
  n.s <- nrow(S)
  
  for (i in 1:n) {
    subtotal <- 0
    x <- X[i, ]
    y.point <- y[i]
    
    for (s in 1:n.s) {
      subtotal <- subtotal + point.likelihood(S[s, ], x, y.point)
    }
    
    total <- total + log(subtotal/n.s)
  }
  
  return(total)
}
