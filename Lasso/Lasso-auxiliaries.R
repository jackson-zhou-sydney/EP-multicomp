
# Auxiliary functions and variables for lasso linear regression

mu.kappa <- 0
sigma.2.kappa <- 10000
lambda <- 0.5
sigma <- 0.1

expec.lnig <- function(A, B, C, D, upper = 100) {
  # Expectation of product of log-normal and inverse gamma densities
  p <- function(x) -A*log(x) - (log(x) - B)^2/(2*C) - D/(2*x^2)
  p.max <- optimise(p, c(0, upper), maximum = TRUE)$objective
  q <- function(x) exp(-A*log(x) - (log(x) - B)^2/(2*C) - D/(2*x^2) - p.max)
  integrate(function(x) q(x)/(x^2), 0, Inf, abs.tol = 0)$value/integrate(q, 0, Inf, abs.tol = 0)$value
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
  expec.inv.sigma.2 <- 1
  expec.a <- rep(1, p)
  
  # Main MFVB loop
  for (iteration in 1:maxit) {
    # Storing old values
    mu.beta.old <- mu.beta
    
    # Update q(beta)
    Q.inv <- XTX + lambda^2*diag(expec.a)
    Q <- solve(Q.inv, tol = 1.0E-99)
    mu.beta <- as.vector(Q%*%XTy)
    Sigma.beta <- Q/expec.inv.sigma.2
    
    # Update q(sigma)
    expec.inv.sigma.2 <- expec.lnig(n + p + 1, mu.kappa, sigma.2.kappa, 
                                    sum((y - X%*%mu.beta)^2) + 
                                    lambda^2*sum(expec.a*mu.beta^2) + 
                                    sum(diag(Q.inv%*%Sigma.beta)))
    
    # Update q(a)
    expec.a <- sqrt(1/(expec.inv.sigma.2*(mu.beta^2 + diag(Sigma.beta))))/lambda
    
    # Checking for convergence
    delta <- norm(mu.beta - mu.beta.old, "2")
    if (verbose) print(paste0("Iteration ", iteration, ", current delta is: ", delta))
    if (delta < tol) break
  }
  
  # Returning only beta
  return(list(mu = mu.beta, Sigma = Sigma.beta))
}

I.r.1 <- function(y, m, V, eta, mult, maxEval, tol) {
  # Hybrid integrals for likelihood sites
  m.1 <- m[1]
  m.2 <- m[2]
  
  V.inv <- solve(V)
  V.11 <- V.inv[1, 1]
  V.12 <- V.inv[1, 2]
  V.22 <- V.inv[2, 2]
  
  abc <- function(x) {
    c(V.11 + eta/exp(2*x),
      2*(V.12*(x - m.2) - V.11*m.1) - eta*2*y/exp(2*x),
      V.11*m.1^2 + 2*V.12*m.1*(m.2 - x) + V.22*(x - m.2)^2 + eta*(2*x + y^2/exp(2*x)))
  }
  
  lb <- m.2 - mult*sqrt(V[2, 2])
  ub <- m.2 + mult*sqrt(V[2, 2])
  
  ret.0 <- integrate(Vectorize(function(x) GI.0(abc(x))), lb, ub)$value
  ret.11 <- integrate(Vectorize(function(x) GI.1(abc(x))), lb, ub)$value
  ret.12 <- integrate(Vectorize(function(x) x*GI.0(abc(x))), lb, ub)$value
  ret.211 <- integrate(Vectorize(function(x) GI.2(abc(x))), lb, ub)$value
  ret.212 <- integrate(Vectorize(function(x) x*GI.1(abc(x))), lb, ub)$value
  ret.222 <- integrate(Vectorize(function(x) x^2*GI.0(abc(x))), lb, ub)$value
  
  list(I.0 = ret.0/sqrt(det(2*pi*V)), 
       I.1 = c(ret.11, ret.12)/sqrt(det(2*pi*V)), 
       I.2 = matrix(c(ret.211, ret.212, ret.212, ret.222)/sqrt(det(2*pi*V)), nrow = 2))
}

I.r.2 <- function(lambda, m, V, eta, mult, maxEval, tol) {
  # Hybrid integrals for Laplace-based prior sites
  m.1 <- m[1]
  m.2 <- m[2]
  
  V.inv <- solve(V)
  V.11 <- V.inv[1, 1]
  V.12 <- V.inv[1, 2]
  V.22 <- V.inv[2, 2]
  
  abc <- function(x, lower) {
    c(V.11,
      2*(V.12*(x - m.2) - V.11*m.1) + eta*ifelse(lower, -1, 1)*2*lambda/exp(x),
      V.11*m.1^2 + 2*V.12*m.1*(m.2 - x) + V.22*(x - m.2)^2 + eta*2*x)
  }
  
  lb <- m.2 - mult*sqrt(V[2, 2])
  ub <- m.2 + mult*sqrt(V[2, 2])
  
  ret.0 <- integrate(Vectorize(function(x) TGI.minus.0(abc(x, T), 0) + TGI.plus.0(abc(x, F), 0)), lb, ub)$value
  ret.11 <- integrate(Vectorize(function(x) TGI.minus.1(abc(x, T), 0) + TGI.plus.1(abc(x, F), 0)), lb, ub)$value
  ret.12 <- integrate(Vectorize(function(x) x*(TGI.minus.0(abc(x, T), 0) + TGI.plus.0(abc(x, F), 0))), lb, ub)$value
  ret.211 <- integrate(Vectorize(function(x) TGI.minus.2(abc(x, T), 0) + TGI.plus.2(abc(x, F), 0)), lb, ub)$value
  ret.212 <- integrate(Vectorize(function(x) x*(TGI.minus.1(abc(x, T), 0) + TGI.plus.1(abc(x, F), 0))), lb, ub)$value
  ret.222 <- integrate(Vectorize(function(x) x^2*(TGI.minus.0(abc(x, T), 0) + TGI.plus.0(abc(x, F), 0))), lb, ub)$value
  
  list(I.0 = ret.0/sqrt(det(2*pi*V)), 
       I.1 = c(ret.11, ret.12)/sqrt(det(2*pi*V)), 
       I.2 = matrix(c(ret.211, ret.212, ret.212, ret.222)/sqrt(det(2*pi*V)), nrow = 2))
}

ep.approx <- function(X, y, mu.kappa, sigma.2.kappa, 
                      lambda, eta, alpha, Q.star, r.star, prec, 
                      min.passes, max.passes, tol.factor, stop.factor, 
                      abs.thresh, rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian lasso linear regression
  n <- nrow(X)
  p <- ncol(X)
  stop.ep <- F
  
  # Parameter initialisation
  Q.star.values <- array(dim = c(2, 2, n + p))
  r.star.values <- matrix(nrow = n + p, ncol = 2)
  
  Q.p <- matrix(0, nrow = p + 1, ncol = p + 1)
  Q.p[p + 1, p + 1] <- 1/sigma.2.kappa
  r.p <- numeric(p + 1)
  r.p[p + 1] <- mu.kappa/sigma.2.kappa
  
  Q.sum <- prec*diag(p + 1) + Q.p
  r.sum <- r.p
  
  for (i in 1:(n + p)) {
    W <- matrix(0, nrow = p + 1, ncol = 2)
    if (i <= n) {
      W[1:p, 1] <- X[i, ]
      W[p + 1, 2] <- 1
    } else {
      W[i - n, 1] <- 1
      W[p + 1, 2] <- 1
    }
    
    Q.star.values[, , i] <- Q.star
    r.star.values[i, ] <- r.star
    Q.sum <- Q.sum + force.sym(W%*%Q.star.values[, , i]%*%t(W))
    r.sum <- r.sum + W%*%r.star.values[i, ]
  }
  
  # Delta initialisation
  deltas <- matrix(-1, nrow = max.passes*(n + p), ncol = 5)
  colnames(deltas) <- c("index", "iteration", "i", "delta", "skip")
  prev.max.delta <- Inf
  prev.med.delta <- Inf
  pcount <- 0
  index <- 1
  
  # Main EP loop
  for (iteration in 1:max.passes) {
    
    if (verbose) print(paste0("---- Current iteration: ", iteration, " ----"))
    
    for (i in sample(1:(n + p))) {
      # Setting up
      W <- matrix(0, nrow = p + 1, ncol = 2)
      if (i <= n) {
        W[1:p, 1] <- X[i, ]
        W[p + 1, 2] <- 1
      } else {
        W[i - n, 1] <- 1
        W[p + 1, 2] <- 1
      }
      
      Q.cavity <- Q.sum - eta*force.sym(W%*%Q.star.values[, , i]%*%t(W))
      Q.cavity.inv <- tryCatch(force.sym(solve(Q.cavity)), error = err)
      if (!is.matrix(Q.cavity.inv)) {stop.ep <- T; break}
      r.cavity <- r.sum - eta*W%*%r.star.values[i, ]
      mu.cavity <- Q.cavity.inv%*%r.cavity
      Sigma.cavity <- Q.cavity.inv
      m <- t(W)%*%mu.cavity
      V <- force.sym(t(W)%*%Sigma.cavity%*%W)
      if (det(V) < 0) {
        print(paste0("Warning: bad V at i = ", i))
        deltas[index, ] <- c(index, iteration, i, NA, 1)
        index <- index + 1
        next
      }
      U <- Sigma.cavity%*%W
      
      # Computing function values, gradients and Hessians at 0
      if (i <= n) {
        I.r.res <- tryCatch(I.r.1(y[i], m, V, eta, mult = 5, maxEval = 0, tol = 0.0001), error = err)
      } else {
        I.r.res <- tryCatch(I.r.2(lambda, m, V, eta, mult = 5, maxEval = 0, tol = 0.0001), error = err)
      }
      
      if (!is.list(I.r.res)) {
        print(paste0("Warning: error in hybrid integral at i = ", i))
        deltas[index, ] <- c(index, iteration, i, NA, 1)
        index <- index + 1
        next
      }
      
      I.0 <- I.r.res$I.0
      I.1 <- I.r.res$I.1
      I.2 <- I.r.res$I.2
      
      f.0 <- 1
      f.grad.0 <- mu.cavity
      f.hess.0 <- mu.cavity%*%t(mu.cavity) + Sigma.cavity
      
      g.0 <- I.0
      g.grad.0 <- U%*%solve(V)%*%(I.1 - m*I.0)
      g.hess.0 <- U%*%solve(V)%*%(I.2 - I.1%*%t(m) - m%*%t(I.1) + m%*%t(m)*I.0)%*%solve(V)%*%t(U) - I.0*U%*%solve(V)%*%t(U)
      
      # Computing hybrid moments
      z.hybrid <- f.0*g.0
      mu.hybrid <- (f.0*g.grad.0 + f.grad.0*g.0)/z.hybrid
      Sigma.hybrid <- (f.hess.0*g.0 + f.grad.0%*%t(g.grad.0) + g.grad.0%*%t(f.grad.0) + f.0*g.hess.0)/z.hybrid - mu.hybrid%*%t(mu.hybrid)
      Sigma.hybrid.inv <- tryCatch(force.sym(solve(Sigma.hybrid)), error = err)
      if (!is.matrix(Sigma.hybrid.inv)) {stop.ep <- T; break}
      
      # Moment matching and calculating deltas
      Q.updated <- (Sigma.hybrid.inv - Q.cavity)/eta
      r.updated <- (Sigma.hybrid.inv%*%mu.hybrid - r.cavity)/eta
      
      W.r <- rowSums(W)
      Q.ratio <- Q.updated/(W.r%*%t(W.r))
      r.ratio <- r.updated/W.r
      
      if (i <= n) {
        Q.star.updated <- force.sym(block.mean(Q.ratio, list(1:p, p + 1)))
        r.star.updated <- c(mean(r.ratio[1:p]), r.ratio[p + 1])
      } else {
        Q.star.updated <- force.sym(Q.ratio[c(i - n, p + 1), c(i - n, p + 1)])
        r.star.updated <- r.ratio[c(i - n, p + 1)]
      }
      
      delta <- max(norm(r.star.updated - r.star.values[i, ], "2"), norm(Q.star.updated - Q.star.values[, , i], "F"))
      
      if (is.na(delta) || delta > tol.factor*prev.med.delta || delta > delta.limit) {
        deltas[index, ] <- c(index, iteration, i, delta, 1)
        index <- index + 1
        next
      } else {
        deltas[index, ] <- c(index, iteration, i, delta, 0)
        index <- index + 1
      }
      
      Q.star.new <- (1 - alpha)*Q.star.values[, , i] + alpha*Q.star.updated
      r.star.new <- (1 - alpha)*r.star.values[i, ] + alpha*r.star.updated
      Q.sum <- Q.sum - force.sym(W%*%Q.star.values[, , i]%*%t(W)) + force.sym(W%*%Q.star.new%*%t(W))
      r.sum <- r.sum - W%*%r.star.values[i, ] + W%*%r.star.new
      Q.star.values[, , i] <- Q.star.new
      r.star.values[i, ] <- r.star.new
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
  Q.final <- Q.sum
  Q.final.inv <- solve(Q.final)
  r.final <- r.sum
  return(list(mu = Q.final.inv%*%r.final, Sigma = Q.final.inv, deltas = deltas[deltas[, "index"] != -1, ]))
}
