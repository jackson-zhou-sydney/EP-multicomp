
# Auxiliary functions and variables for heteroscedastic linear regression

sourceCpp("Hetero/Methods/EP.cpp")
sourceCpp("Hetero/Methods/EP-2D.cpp")
source("Hetero/Methods/GVB.R")

sigma.2.beta.1 <- 1
sigma.2.beta.2 <- 0.01

sim.settings <- list(c(n = 200, p.1 = 40, p.2 = 10),
                     c(n = 200, p.1 = 20, p.2 = 20),
                     c(n = 200, p.1 = 10, p.2 = 40))

sim.labels <- c("1" = "n = 200, p_1 = 40, p_2 = 10",
                "2" = "n = 200, p_1 = 20, p_2 = 20",
                "3" = "n = 200, p_1 = 10, p_2 = 40")

bench.settings <- list(c(n = 40, p.1 = 2, p.2 = 2),
                       c(n = 725, p.1 = 6, p.2 = 2),
                       c(n = 125, p.1 = 5, p.2 = 5),
                       c(n = 17268, p.1 = 172, p.2 = 172))

bench.labels <- c("1" = "Food (n = 40, p_1 = 2, p_2 = 2)",
                  "2" = "Salary (n = 725, p_1 = 6, p_2 = 2)",
                  "3" = "Sniffer (n = 125, p_1 = 5, p_2 = 5)",
                  "4" = "Energy (n = 17268, p_1 = 172, p_2 = 172)")

point.likelihood <- function(theta, x.1, x.2, y) {
  # Likelihood evaluated at a point
  p.1 <- length(x.1)
  p.2 <- length(x.2)
  
  beta.1 <- theta[1:p.1]
  beta.2 <- theta[(p.1 + 1):length(theta)]
  
  as.numeric(exp(-t(x.2)%*%beta.2 - (y - t(x.1)%*%beta.1)^2/as.numeric(2*exp(2*t(x.2)%*%beta.2)))/sqrt(2*pi))
}

lppd <- function(X.1, X.2, y, S) {
  # Compute the log pointwise predictive density
  # S is a matrix of posterior samples
  total <- 0
  n <- nrow(X.1)
  n.s <- nrow(S)
  
  for (i in 1:n) {
    subtotal <- 0
    x.1 <- X.1[i, ]
    x.2 <- X.2[i, ]
    y.point <- y[i]
    
    for (s in 1:n.s) {
      subtotal <- subtotal + point.likelihood(S[s, ], x.1, x.2, y.point)
    }
    
    total <- total + log(subtotal/n.s)
  }
  
  return(total)
}
