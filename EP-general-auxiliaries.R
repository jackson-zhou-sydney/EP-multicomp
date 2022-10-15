
# General auxiliary functions and variables for expectation propagation

err <- function(e) {
  # Return NA on error
  return(NA)
}

force.sym <- function(m) {
  # Force matrix to be symmetric
  m.sym <- m
  m.sym[lower.tri(m.sym)] <- t(m.sym)[lower.tri(m.sym)]
  return(m.sym)
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

TGI.lower.0 <- function(x, y) {
  # Truncated Gaussian integral (lower, 0th raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  exp(log(sqrt(2*pi/a)) - 0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, log.p = T))
}

TGI.lower.1 <- function(x, y) {
  # Truncated Gaussian integral (lower, 1st raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*(-(b/(2*a))*exp(-0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, log.p = T)) - (1/sqrt(a))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))
}

TGI.lower.2 <- function(x, y) {
  # Truncated Gaussian integral (lower, 2nd raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*((b^2/(4*a) + 1)*exp(-0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, log.p = T)) + (b/(2*sqrt(a)) - sqrt(a)*y)*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))/a
}

TGI.upper.0 <- function(x, y) {
  # Truncated Gaussian integral (upper, 0th raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  exp(log(sqrt(2*pi/a)) - 0.5*(c - b^2/(4*a)) + pnorm(b/(2*sqrt(a)) + sqrt(a)*y, lower.tail = F, log.p = T))
}

TGI.upper.1 <- function(x, y) {
  # Truncated Gaussian integral (upper, 1st raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*(-(b/(2*a))*exp(-0.5*(c - b^2/(4*a)) + pnorm(-b/(2*sqrt(a)) - sqrt(a)*y, log.p = T)) + (1/sqrt(a))*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))
}

TGI.upper.2 <- function(x, y) {
  # Truncated Gaussian integral (upper, 2nd raw moment)
  a <- x[1]; b <- x[2]; c <- x[3]
  sqrt(2*pi/a)*((b^2/(4*a) + 1)*exp(-0.5*(c - b^2/(4*a)) + pnorm(-b/(2*sqrt(a)) - sqrt(a)*y, log.p = T)) - (b/(2*sqrt(a)) - sqrt(a)*y)*exp(-0.5*(c - b^2/(4*a)) + dnorm(b/(2*sqrt(a)) + sqrt(a)*y, log = T)))/a
}
