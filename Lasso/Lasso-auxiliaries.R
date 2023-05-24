
# Auxiliary functions and variables for lasso linear regression

sourceCpp("Lasso/EP-approx.cpp")
sourceCpp("Lasso/MFVB-approx.cpp")

mu.kappa <- 0
sigma.2.kappa <- 10000
lambda <- 0.5
kappa <- -1

sim.settings <- list(c(n = 200, p = 40),
                     c(n = 40, p = 40),
                     c(n = 10, p = 40))

bench.settings <- list(c(n = 442, p = 11),
                       c(n = 97, p = 9),
                       c(n = 120, p = 201))

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

expec.lnig <- function(A, B, C, D, fun, radius = 100) {
  # Expectation of product of reparameterised log-normal and inverse gamma densities
  p <- function(x) -A*x - (x - B)^2/(2*C) - D/(2*exp(2*x))
  p.max <- optimise(p, c(-radius, radius), maximum = TRUE)$objective
  q <- function(x) exp(-A*x - (x - B)^2/(2*C) - D/(2*exp(2*x)) - p.max)
  
  if (fun == "x") {
    r <- function(x) exp(-A*x - (x - B)^2/(2*C) - D/(2*exp(2*x)) - p.max)*x
  } else if (fun == "x^2") {
    r <- function(x) exp(-A*x - (x - B)^2/(2*C) - D/(2*exp(2*x)) - p.max)*x^2
  } else if (fun == "1/exp(2*x)") {
    r <- function(x) exp(-A*x - (x - B)^2/(2*C) - D/(2*exp(2*x)) - p.max - 2*x)
  } else {
    stop("fun must be one of: x, x^2, or 1/exp(2*x)")
  }
  
  integrate(r, -Inf, Inf, abs.tol = 0)$value/integrate(q, -Inf, Inf, abs.tol = 0)$value
}

mfvb.approx <- function(X, y, mu.kappa, sigma.2.kappa,
                        lambda, maxit, tol, verbose) {
  # MFVB for Bayesian lasso linear regression
  n <- nrow(X)
  p <- ncol(X)
  
  XTX <- t(X)%*%X
  XTy <- t(X)%*%y
  
  # Initialisations
  mu.beta <- rep(0, p)
  Sigma.beta <- diag(p)
  expec.ie.2.kappa <- 1
  expec.a <- rep(1, p)
  
  # Main MFVB loop
  for (iteration in 1:maxit) {
    # Storing old values
    mu.beta.old <- mu.beta
    
    # Update q(beta)
    Q.inv <- XTX + lambda^2*diag(expec.a)
    Q <- solve(Q.inv, tol = 1.0E-99)
    mu.beta <- as.vector(Q%*%XTy)
    Sigma.beta <- Q/expec.ie.2.kappa
    
    # Update q(kappa)
    expec.ie.2.kappa <- expec.lnig(n + p, mu.kappa, sigma.2.kappa,
                                   sum((y - X%*%mu.beta)^2) + 
                                   lambda^2*sum(expec.a*mu.beta^2) + 
                                   sum(diag(Q.inv%*%Sigma.beta)),
                                   fun = "1/exp(2*x)")
    
    # Update q(a)
    expec.a <- sqrt(1/(expec.ie.2.kappa*(mu.beta^2 + diag(Sigma.beta))))/lambda
    
    # Checking for convergence
    delta <- norm(mu.beta - mu.beta.old, "2")
    if (verbose) print(paste0("Iteration ", iteration, ", current delta is: ", delta))
    if (delta < tol) break
  }
  
  # Returning parameters
  mu.kappa.q <- expec.lnig(n + p, mu.kappa, sigma.2.kappa,
                           sum((y - X%*%mu.beta)^2) + 
                           lambda^2*sum(expec.a*mu.beta^2) + 
                           sum(diag(Q.inv%*%Sigma.beta)),
                           fun = "x")
  
  sigma.2.kappa.q <- expec.lnig(n + p, mu.kappa, sigma.2.kappa,
                                sum((y - X%*%mu.beta)^2) + 
                                lambda^2*sum(expec.a*mu.beta^2) + 
                                sum(diag(Q.inv%*%Sigma.beta)),
                                fun = "x^2") - mu.kappa.q^2
  
  mu.theta <- numeric(p + 1)
  mu.theta[1:p] <- mu.beta
  mu.theta[p + 1] <- mu.kappa.q
  
  Sigma.theta <- matrix(0, nrow = p + 1, ncol = p + 1)
  Sigma.theta[1:p, 1:p] <- Sigma.beta
  Sigma.theta[p + 1, p + 1] <- sigma.2.kappa.q
  
  return(list(mu = mu.theta, Sigma = Sigma.theta))
}

h.mom.1 <- function(y, mu, Sigma, eta, mult, lb.min, ub.max, length) {
  # Hybrid moments for likelihood sites
  mu.1 <- mu[1]
  mu.2 <- mu[2]
  
  Q <- sym(solve(Sigma))
  Q.11 <- Q[1, 1]
  Q.12 <- Q[1, 2]
  Q.22 <- Q[2, 2]
  
  abc <- function(x) {
    c(Q.11 + eta/exp(2*x),
      2*(Q.12*(x - mu.2) - Q.11*mu.1) - eta*2*y/exp(2*x),
      Q.11*mu.1^2 + 2*Q.12*mu.1*(mu.2 - x) + Q.22*(x - mu.2)^2 + eta*(2*x + y^2/exp(2*x) - 2*mu.2 - (y - mu.1)^2/exp(2*mu.2)))
  }
  
  lb <- max(mu.2 - mult*sqrt(Sigma[2, 2]), lb.min)
  ub <- min(mu.2 + mult*sqrt(Sigma[2, 2]), ub.max)
  
  x.values <- seq(from = lb, to = ub, length = length)
  y.matrix <- matrix(nrow = 6, ncol = length)
  
  for (i in 1:length(x.values)) {
    x <- x.values[i]
    abc.value <- abc(x)
    y.matrix[1, i] <- GI.0(abc.value)
    y.matrix[2, i] <- GI.1(abc.value)
    y.matrix[3, i] <- x*y.matrix[1, i]
    y.matrix[4, i] <- GI.2(abc.value)
    y.matrix[5, i] <- x*y.matrix[2, i]
    y.matrix[6, i] <- x^2*y.matrix[1, i]
  }
  
  ret.0 <- trapz(x.values, y.matrix[1, ])
  ret.11 <- trapz(x.values, y.matrix[2, ])
  ret.12 <- trapz(x.values, y.matrix[3, ])
  ret.211 <- trapz(x.values, y.matrix[4, ])
  ret.212 <- trapz(x.values, y.matrix[5, ])
  ret.222 <- trapz(x.values, y.matrix[6, ])
  
  mu.h <- c(ret.11, ret.12)/ret.0
  Sigma.h <- sym(matrix(c(ret.211, ret.212, ret.212, ret.222)/ret.0, nrow = 2) - mu.h%*%t(mu.h))
  return(list(mu.h = mu.h, Sigma.h = Sigma.h))
}

h.mom.2 <- function(lambda, mu, Sigma, eta, mult, lb.min, ub.max, length) {
  # Hybrid integrals for Laplace-based prior sites
  mu.1 <- mu[1]
  mu.2 <- mu[2]
  
  Q <- sym(solve(Sigma))
  Q.11 <- Q[1, 1]
  Q.12 <- Q[1, 2]
  Q.22 <- Q[2, 2]
  
  abc <- function(x, lower) {
    c(Q.11,
      2*(Q.12*(x - mu.2) - Q.11*mu.1) + eta*ifelse(lower, -1, 1)*2*lambda/exp(x),
      Q.11*mu.1^2 + 2*Q.12*mu.1*(mu.2 - x) + Q.22*(x - mu.2)^2 + eta*2*x)
  }
  
  lb <- max(mu.2 - mult*sqrt(Sigma[2, 2]), lb.min)
  ub <- min(mu.2 + mult*sqrt(Sigma[2, 2]), ub.max)
  
  x.values <- seq(from = lb, to = ub, length = length)
  y.matrix <- matrix(nrow = 6, ncol = length)
  
  for (i in 1:length(x.values)) {
    x <- x.values[i]
    abc.l <- abc(x, T)
    abc.u <- abc(x, F)
    y.matrix[1, i] <- TGI.minus.0(abc.l, 0) + TGI.plus.0(abc.u, 0)
    y.matrix[2, i] <- TGI.minus.1(abc.l, 0) + TGI.plus.1(abc.u, 0)
    y.matrix[3, i] <- x*y.matrix[1, i]
    y.matrix[4, i] <- TGI.minus.2(abc.l, 0) + TGI.plus.2(abc.u, 0)
    y.matrix[5, i] <- x*y.matrix[2, i]
    y.matrix[6, i] <- x^2*y.matrix[1, i]
  }
  
  ret.0 <- trapz(x.values, y.matrix[1, ])
  ret.11 <- trapz(x.values, y.matrix[2, ])
  ret.12 <- trapz(x.values, y.matrix[3, ])
  ret.211 <- trapz(x.values, y.matrix[4, ])
  ret.212 <- trapz(x.values, y.matrix[5, ])
  ret.222 <- trapz(x.values, y.matrix[6, ])
  
  mu.h <- c(ret.11, ret.12)/ret.0
  Sigma.h <- sym(matrix(c(ret.211, ret.212, ret.212, ret.222)/ret.0, nrow = 2) - mu.h%*%t(mu.h))
  return(list(mu.h = mu.h, Sigma.h = Sigma.h))
}

ep.approx <- function(X, y, sigma.2.kappa, mu.kappa, lambda, eta, alpha, Q.star.init, r.star.init,
                      min.passes, max.passes, thresh, verbose) {
  # Dampened power EP for Bayesian lasso linear regression
  n <- nrow(X)
  p <- ncol(X)
  
  # Parameter initialisation
  Q.star.values <- array(dim = c(2, 2, n + p))
  r.star.values <- matrix(nrow = n + p, ncol = 2)
  
  Q.p <- matrix(0, nrow = p + 1, ncol = p + 1)
  Q.p[p + 1, p + 1] <- 1/sigma.2.kappa
  r.p <- numeric(p + 1)
  r.p[p + 1] <- mu.kappa/sigma.2.kappa
  
  Q.dot <- Q.p
  r.dot <- r.p
  
  for (k in 1:(n + p)) {
    A <- matrix(0, nrow = p + 1, ncol = 2)
    if (k <= n) {
      A[1:p, 1] <- X[k, ]
      A[p + 1, 2] <- 1
    } else {
      A[k - n, 1] <- 1
      A[p + 1, 2] <- 1
    }
    Q.star.values[, , k] <- sym(Q.star.init)
    r.star.values[k, ] <- r.star.init
    Q.dot <- Q.dot + sym(A%*%Q.star.values[, , k]%*%t(A))
    r.dot <- r.dot + A%*%r.star.values[k, ]
  }
  
  Sigma.dot <- sym(solve(Q.dot))
  mu.dot <- Sigma.dot%*%r.dot
  
  # Main EP loop
  for (pass in 1:max.passes) {
    
    # Delta initialisation
    deltas.Q <- rep(NA, n + p)
    deltas.r <- rep(NA, n + p)
    
    if (verbose) print(paste0("---- Current pass: ", pass, " ----"))
    
    for (k in sample(1:(n + p))) {
      
      A <- matrix(0, nrow = p + 1, ncol = 2)
      if (k <= n) {
        A[1:p, 1] <- X[k, ]
        A[p + 1, 2] <- 1
      } else {
        A[k - n, 1] <- 1
        A[p + 1, 2] <- 1
      }
      Sigma.dot.star <- sym(t(A)%*%Sigma.dot%*%A); mu.dot.star <- t(A)%*%mu.dot
      Q.dot.star <- sym(solve(Sigma.dot.star)); r.dot.star <- Q.dot.star%*%mu.dot.star
      Q.m.star <- Q.dot.star - eta*Q.star.values[, , k]; r.m.star <- r.dot.star - eta*r.star.values[k, ]
      Sigma.m.star <- sym(solve(Q.m.star)); mu.m.star <- Sigma.m.star%*%r.m.star
      
      if (k <= n) {
        h.mom.res <- tryCatch(h.mom.1(y[k], mu.m.star, Sigma.m.star, eta, mult = 5, lb.min = -5, ub.max = Inf, length = 50), error = err)
      } else {
        h.mom.res <- tryCatch(h.mom.2(lambda, mu.m.star, Sigma.m.star, eta, mult = 5, lb.min = -10, ub.max = Inf, length = 50), error = err)
      }
      if (!is.list(h.mom.res)) {
        print(paste0("Warning: error in hybrid integral at k = ", k))
        next
      }
      Sigma.h.star <- h.mom.res$Sigma.h; mu.h.star <- h.mom.res$mu.h
      Q.h.star <- sym(solve(Sigma.h.star)); r.h.star <- Q.h.star%*%mu.h.star
      Q.star.tilde <- (1 - alpha)*Q.star.values[, , k] + (alpha/eta)*(Q.h.star - Q.m.star)
      r.star.tilde <- (1 - alpha)*r.star.values[k, ] + (alpha/eta)*(r.h.star - r.m.star)
      Q.star.tilde.d <- Q.star.tilde - Q.star.values[, , k]; r.star.tilde.d <- r.star.tilde - r.star.values[k, ]
      
      deltas.Q[k] <- norm(Q.star.tilde.d, "F")
      deltas.r[k] <- norm(r.star.tilde.d, "2")
      
      Q.dot <- Q.dot + sym(A%*%Q.star.tilde.d%*%t(A)); r.dot <- r.dot + A%*%(r.star.tilde.d)
      Sigma.dot <- Sigma.dot - sym((Sigma.dot%*%A)%*%solve(solve(Q.star.tilde.d) + Sigma.dot.star)%*%t(Sigma.dot%*%A)); mu.dot <- Sigma.dot%*%r.dot
      Q.star.values[, , k] <- Q.star.tilde; r.star.values[k, ] <- r.star.tilde
      
    }
    
    if (pass == 1) {
      # Base maximum deltas
      bmd.Q <- max(deltas.Q, na.rm = T)
      bmd.r <- max(deltas.r, na.rm = T)
      
      if (verbose) {
        print(paste0("Maximum delta for Q: ", round(bmd.Q, 2)))
        print(paste0("Maximum delta for r: ", round(bmd.r, 2)))
      }
    } else {
      md.Q <- max(deltas.Q, na.rm = T)
      md.r <- max(deltas.r, na.rm = T)
      
      if (verbose) {
        print(paste0("Maximum delta for Q: ", round(md.Q, 2)))
        print(paste0("Maximum delta for r: ", round(md.r, 2)))
      }
      
      if (md.Q <= thresh*bmd.Q && md.r <= thresh*bmd.r && pass > min.passes) {
        if (verbose) print("EP has converged; stopping EP")
        break
      }
    }
    
  }
  
  # Returning in original parameterisation
  return(list(mu = mu.dot, Sigma = Sigma.dot))
}
