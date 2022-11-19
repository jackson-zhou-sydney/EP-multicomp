
# Auxiliary functions and variables for lasso linear regression

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

I.h.r.1 <- function(y, mu, Sigma, eta, mult, lb.min, ub.max) {
  # Hybrid integrals for likelihood sites
  mu.1 <- mu[1]
  mu.2 <- mu[2]
  
  Q <- solve(Sigma)
  Q.11 <- Q[1, 1]
  Q.12 <- Q[1, 2]
  Q.22 <- Q[2, 2]
  
  abc <- function(x) {
    c(Q.11 + eta/exp(2*x),
      2*(Q.12*(x - mu.2) - Q.11*mu.1) - eta*2*y/exp(2*x),
      Q.11*mu.1^2 + 2*Q.12*mu.1*(mu.2 - x) + Q.22*(x - mu.2)^2 + eta*(2*x + y^2/exp(2*x)))
  }
  
  lb <- max(mu.2 - mult*sqrt(Sigma[2, 2]), lb.min)
  ub <- min(mu.2 + mult*sqrt(Sigma[2, 2]), ub.max)
  
  ret.0 <- integrate(Vectorize(function(x) GI.0(abc(x))), lb, ub, abs.tol = 0)$value
  ret.11 <- integrate(Vectorize(function(x) GI.1(abc(x))), lb, ub, abs.tol = 0)$value
  ret.12 <- integrate(Vectorize(function(x) x*GI.0(abc(x))), lb, ub, abs.tol = 0)$value
  ret.211 <- integrate(Vectorize(function(x) GI.2(abc(x))), lb, ub, abs.tol = 0)$value
  ret.212 <- integrate(Vectorize(function(x) x*GI.1(abc(x))), lb, ub, abs.tol = 0)$value
  ret.222 <- integrate(Vectorize(function(x) x^2*GI.0(abc(x))), lb, ub, abs.tol = 0)$value
  
  list(I.h.0 = ret.0/sqrt(det(2*pi*Sigma)), 
       I.h.1 = c(ret.11, ret.12)/sqrt(det(2*pi*Sigma)), 
       I.h.2 = matrix(c(ret.211, ret.212, ret.212, ret.222)/sqrt(det(2*pi*Sigma)), nrow = 2))
}

I.h.r.2 <- function(lambda, mu, Sigma, eta, mult, lb.min, ub.max) {
  # Hybrid integrals for Laplace-based prior sites
  mu.1 <- mu[1]
  mu.2 <- mu[2]
  
  Q <- solve(Sigma)
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
  
  ret.0 <- integrate(Vectorize(function(x) TGI.minus.0(abc(x, T), 0) + TGI.plus.0(abc(x, F), 0)), lb, ub, abs.tol = 0)$value
  ret.11 <- integrate(Vectorize(function(x) TGI.minus.1(abc(x, T), 0) + TGI.plus.1(abc(x, F), 0)), lb, ub, abs.tol = 0)$value
  ret.12 <- integrate(Vectorize(function(x) x*(TGI.minus.0(abc(x, T), 0) + TGI.plus.0(abc(x, F), 0))), lb, ub, abs.tol = 0)$value
  ret.211 <- integrate(Vectorize(function(x) TGI.minus.2(abc(x, T), 0) + TGI.plus.2(abc(x, F), 0)), lb, ub, abs.tol = 0)$value
  ret.212 <- integrate(Vectorize(function(x) x*(TGI.minus.1(abc(x, T), 0) + TGI.plus.1(abc(x, F), 0))), lb, ub, abs.tol = 0)$value
  ret.222 <- integrate(Vectorize(function(x) x^2*(TGI.minus.0(abc(x, T), 0) + TGI.plus.0(abc(x, F), 0))), lb, ub, abs.tol = 0)$value
  
  list(I.h.0 = ret.0/sqrt(det(2*pi*Sigma)), 
       I.h.1 = c(ret.11, ret.12)/sqrt(det(2*pi*Sigma)), 
       I.h.2 = matrix(c(ret.211, ret.212, ret.212, ret.222)/sqrt(det(2*pi*Sigma)), nrow = 2))
}

ep.approx <- function(X, y, mu.kappa, sigma.2.kappa, 
                      lambda, eta, alpha, Q.star.init, r.star.init, offset,
                      min.passes, max.passes, tol.factor, stop.factor,
                      abs.thresh, rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian lasso linear regression
  n <- nrow(X)
  p <- ncol(X)
  stop.ep <- F
  pcount <- 0
  
  # Parameter initialisation
  Q.star.values <- array(dim = c(2, 2, n + p))
  r.star.values <- matrix(nrow = n + p, ncol = 2)
  
  Q.p <- matrix(0, nrow = p + 1, ncol = p + 1)
  Q.p[p + 1, p + 1] <- 1/sigma.2.kappa
  r.p <- numeric(p + 1)
  r.p[p + 1] <- mu.kappa/sigma.2.kappa
  
  Q.dot <- sym(offset) + Q.p
  r.dot <- r.p
  
  for (i in 1:(n + p)) {
    A <- matrix(0, nrow = p + 1, ncol = 2)
    if (i <= n) {
      A[1:p, 1] <- X[i, ]
      A[p + 1, 2] <- 1
    } else {
      A[i - n, 1] <- 1
      A[p + 1, 2] <- 1
    }
    Q.star.values[, , i] <- sym(Q.star.init)
    r.star.values[i, ] <- r.star.init
    Q.dot <- Q.dot + sym(A%*%Q.star.values[, , i]%*%t(A))
    r.dot <- r.dot + A%*%r.star.values[i, ]
  }
  
  Sigma.dot <- sym(solve(Q.dot))
  mu.dot <- Sigma.dot%*%r.dot
  
  # Delta initialisation
  deltas <- matrix(-1, nrow = max.passes*(n + p), ncol = 5)
  colnames(deltas) <- c("index", "iteration", "i", "delta", "skip")
  prev.max.delta <- Inf
  prev.med.delta <- Inf
  index <- 1
  
  # Main EP loop
  for (iteration in 1:max.passes) {
    
    if (verbose) print(paste0("---- Current iteration: ", iteration, " ----"))
    
    for (i in sample(1:(n + p))) {
      
      A <- matrix(0, nrow = p + 1, ncol = 2)
      if (i <= n) {
        A[1:p, 1] <- X[i, ]
        A[p + 1, 2] <- 1
      } else {
        A[i - n, 1] <- 1
        A[p + 1, 2] <- 1
      }
      Sigma.dot.star <- sym(t(A)%*%Sigma.dot%*%A); mu.dot.star <- t(A)%*%mu.dot
      Q.dot.star <- sym(solve(Sigma.dot.star)); r.dot.star <- Q.dot.star%*%mu.dot.star
      Q.m.star <- Q.dot.star - eta*Q.star.values[, , i]; r.m.star <- r.dot.star - eta*r.star.values[i, ]
      Sigma.m.star <- sym(solve(Q.m.star)); mu.m.star <- Sigma.m.star%*%r.m.star
      
      if (i <= n) {
        I.h.r.res <- tryCatch(I.h.r.1(y[i], mu.m.star, Sigma.m.star, eta, mult = 5, lb.min = -5, ub.max = Inf), error = err)
      } else {
        I.h.r.res <- tryCatch(I.h.r.2(lambda, mu.m.star, Sigma.m.star, eta, mult = 5, lb.min = -10, ub.max = Inf), error = err)
      }
      if (!is.list(I.h.r.res)) {
        print(paste0("Warning: error in hybrid integral at i = ", i))
        deltas[index, ] <- c(index, iteration, i, NA, 1)
        index <- index + 1
        next
      }
      I.h.0 <- I.h.r.res$I.h.0; I.h.1 <- I.h.r.res$I.h.1; I.h.2 <- I.h.r.res$I.h.2
      
      Sigma.h.star <- sym(I.h.2/I.h.0 - I.h.1%*%t(I.h.1)/I.h.0^2); mu.h.star <- I.h.1/I.h.0
      Q.h.star <- sym(solve(Sigma.h.star)); r.h.star <- Q.h.star%*%mu.h.star
      Q.star.tilde <- (1 - alpha)*Q.star.values[, , i] + (alpha/eta)*(Q.h.star - Q.m.star)
      r.star.tilde <- (1 - alpha)*r.star.values[i, ] + (alpha/eta)*(r.h.star - r.m.star)
      Q.star.tilde.d <- Q.star.tilde - Q.star.values[, , i]; r.star.tilde.d <- r.star.tilde - r.star.values[i, ]
      
      delta <- max(norm(Q.star.tilde.d, "F"), norm(r.star.tilde.d, "2"))
      
      if (is.na(delta) || delta > tol.factor*prev.med.delta || delta > delta.limit) {
        deltas[index, ] <- c(index, iteration, i, delta, 1)
        index <- index + 1
        next
      } else {
        deltas[index, ] <- c(index, iteration, i, delta, 0)
        index <- index + 1
      }
      
      Q.dot <- Q.dot + sym(A%*%Q.star.tilde.d%*%t(A)); r.dot <- r.dot + A%*%(r.star.tilde.d)
      Sigma.dot <- Sigma.dot - sym((Sigma.dot%*%A)%*%solve(solve(Q.star.tilde.d) + Sigma.dot.star)%*%t(Sigma.dot%*%A)); mu.dot <- Sigma.dot%*%r.dot
      Q.star.values[, , i] <- Q.star.tilde; r.star.values[i, ] <- r.star.tilde
      
    }
    
    iteration.deltas <- deltas[deltas[, "iteration"] == iteration & deltas[, "skip"] == 0, ][, "delta"]
    max.delta <- max(iteration.deltas)
    med.delta <- median(iteration.deltas)
    if (verbose) {
      print(paste0("Maximum delta: ", round(max.delta, 2)))
      print(paste0("Median delta: ", round(med.delta, 2)))
    }
    
    if (stop.ep) {print("Too many numerical errors; stopping EP"); break}
    if (max.delta < abs.thresh && iteration > min.passes) {if (verbose) print("EP has converged; stopping EP"); break}
    if (max.delta > stop.factor*prev.max.delta) {print("Unstable deltas; stopping EP"); break}
    if (max.delta > rel.thresh*prev.max.delta) pcount <- pcount + 1 else pcount <- 0
    if (pcount == patience) {print("Out of patience; stopping EP"); break}
    
    prev.max.delta <- max.delta
    prev.med.delta <- med.delta
    
  }
  
  # Returning in original parameterisation
  return(list(mu = mu.dot, Sigma = Sigma.dot, deltas = deltas[deltas[, "index"] != -1, ]))
}
