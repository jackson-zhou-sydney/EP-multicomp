
# Auxiliary functions and variables for quantile regression

sourceCpp("Quantile/Methods/EP.cpp")
sourceCpp("Quantile/Methods/EP-2D.cpp")
sourceCpp("Quantile/Methods/MFVB.cpp")

sigma.2.beta <- 1
sigma.2.kappa <- 0.01
tau <- 0.5
sigma <- 0.2
n.binom <- 10

sim.settings <- list(c(n = 200, p = 40, dist = "normal"),
                     c(n = 200, p = 40, dist = "poisson"),
                     c(n = 200, p = 40, dist = "binomial"))

sim.labels <- c("1" = "Normal\nn = 200, p = 40",
                "2" = "Poisson\nn = 200, p = 40",
                "3" = "Binomial\nn = 200, p = 40")

bench.settings <- list(c(n = 298, p = 2),
                       c(n = 235, p = 2),
                       c(n = 21, p = 4),
                       c(n = 17268, p = 172))

bench.labels <- c("1" = "IgG\nn = 298, p = 2",
                  "2" = "Engel\nn = 235, p = 2",
                  "3" = "Stack\nn = 21, p = 4",
                  "4" = "Energy\nn = 17268, p = 172")

rho <- function(x, tau) {
  # Quantile loss function
  0.5*(abs(x) + (2*tau - 1)*x)
}

point.likelihood <- function(theta, x, y, tau) {
  # Likelihood evaluated at a point
  p <- length(x)
  
  beta <- theta[1:p]
  kappa <- theta[p + 1]
  
  as.numeric(tau*(1 - tau)*exp(-kappa - rho(y - t(x)%*%beta, tau)/exp(kappa)))
}

lppd <- function(X, y, tau, S) {
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
      subtotal <- subtotal + point.likelihood(S[s, ], x, y.point, tau)
    }
    
    total <- total + log(subtotal/n.s)
  }
  
  return(total)
}
