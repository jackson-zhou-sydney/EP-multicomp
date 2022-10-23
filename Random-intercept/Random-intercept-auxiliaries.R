
# Auxiliary functions and variables for random intercept linear regression

log.joint.likelihood <- function(theta, X, Z, y, mu.bk, Sigma.bk) {
  # Log joint likelihood
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  beta <- theta[1:p]
  u <- theta[(p + 1):(p + q)]
  kappa.y <- theta[p + q + 1]
  kappa.u <- theta[p + q + 2]
  bk <- c(beta, kappa.y, kappa.u)
  
  as.numeric(-0.5*(2*n*kappa.y + sum((y - X%*%beta - Z%*%u)^2)/exp(2*kappa.y)) - 0.5*(2*q*kappa.u + sum(u^2)/exp(2*kappa.u)) - 0.5*t(bk - mu.bk)%*%solve(Sigma.bk)%*%(bk - mu.bk))
}

laplace.approx <- function(X, Z, y, mu.bk, Sigma.bk, maxit) {
  # Laplace approximation
  p <- ncol(X)
  q <- ncol(Z)
  
  combined.df <- as.data.frame(cbind(X[, -1], group = as.vector(Z%*%1:q), y))
  combined.df$group <- as.factor(combined.df$group)
  lmer.res <- lmer(eval(parse(text = paste0("y ~ (1 | group) + ", paste0("V", 1:(p - 1), collapse = " + ")))), 
                   data = combined.df)
  start <- c(lmer.res@beta, lmer.res@u, log(as.data.frame(VarCorr(lmer.res))$sdcor[2:1]))
  optim.res <- optim(start, fn = log.joint.likelihood, 
                     method = "BFGS", 
                     control = list(fnscale = -1, maxit = maxit),
                     hessian = T,
                     X = X, Z = Z, y = y, mu.bk = mu.bk, Sigma.bk = Sigma.bk)
  return(list(mu = optim.res$par, Sigma = solve(-optim.res$hessian)))
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

I.r.2 <- function(m, V, eta, mult, maxEval, tol) {
  # Hybrid integrals for random effect sites
  m.1 <- m[1]
  m.2 <- m[2]
  
  V.inv <- solve(V)
  V.11 <- V.inv[1, 1]
  V.12 <- V.inv[1, 2]
  V.22 <- V.inv[2, 2]
  
  abc <- function(x) {
    c(V.11 + eta/exp(2*x),
      2*(V.12*(x - m.2) - V.11*m.1),
      V.11*m.1^2 + 2*V.12*m.1*(m.2 - x) + V.22*(x - m.2)^2 + eta*2*x)
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

ep.approx <- function(X, Z, y, mu.bk, Sigma.bk,
                      eta, alpha, Q.star, r.star, prec,
                      min.passes, max.passes, tol.factor, stop.factor, abs.thresh, 
                      rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian random intercept linear regression
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  stop.ep <- F
  
  # Parameter initialisation
  Q.star.values <- array(dim = c(2, 2, n + q))
  r.star.values <- matrix(nrow = n + q, ncol = 2)
  
  Q.bk <- force.sym(solve(Sigma.bk))
  r.bk <- Q.bk%*%mu.bk
  
  Q.p <- matrix(0, nrow = p + q + 2, ncol = p + q + 2)
  Q.p[1:p, 1:p] <- Q.bk[1:p, 1:p]
  Q.p[1:p, (p + q + 1):(p + q + 2)] <- Q.bk[1:p, (p + 1):(p + 2)]
  Q.p[(p + q + 1):(p + q + 2), 1:p] <- Q.bk[(p + 1):(p + 2), 1:p]
  Q.p[(p + q + 1):(p + q + 2), (p + q + 1):(p + q + 2)] <- Q.bk[(p + 1):(p + 2), (p + 1):(p + 2)]
  r.p <- c(r.bk[1:p], rep(0, q), r.bk[(p + 1):(p + 2)])
  
  Q.sum <- prec*diag(p + q + 2) + Q.p
  r.sum <- r.p
  
  for (i in 1:(n + q)) {
    W <- matrix(0, nrow = p + q + 2, ncol = 2)
    if (i <= n) {
      W[1:p, 1] <- X[i, ]
      W[(p + 1):(p + q), 1] <- Z[i, ]
      W[p + q + 1, 2] <- 1
    } else {
      W[p + i - n, 1] <- 1
      W[p + q + 2, 2] <- 1
    }
    
    Q.star.values[, , i] <- Q.star
    r.star.values[i, ] <- r.star
    Q.sum <- Q.sum + force.sym(W%*%Q.star.values[, , i]%*%t(W))
    r.sum <- r.sum + W%*%r.star.values[i, ]
  }
  
  # Delta initialisation
  deltas <- matrix(-1, nrow = max.passes*(n + q), ncol = 5)
  colnames(deltas) <- c("index", "iteration", "i", "delta", "skip")
  prev.max.delta <- Inf
  prev.med.delta <- Inf
  pcount <- 0
  index <- 1
  
  # Main EP loop
  for (iteration in 1:max.passes) {
    
    if (verbose) print(paste0("---- Current iteration: ", iteration, " ----"))
    
    for (i in sample(1:(n + q))) {
      # Setting up
      W <- matrix(0, nrow = p + q + 2, ncol = 2)
      if (i <= n) {
        W[1:p, 1] <- X[i, ]
        W[(p + 1):(p + q), 1] <- Z[i, ]
        W[p + q + 1, 2] <- 1
      } else {
        W[p + i - n, 1] <- 1
        W[p + q + 2, 2] <- 1
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
        I.r.res <- tryCatch(I.r.2(m, V, eta, mult = 5, maxEval = 0, tol = 0.0001), error = err)
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
      Q.ratio[!is.finite(Q.ratio)] <- 0 # Ensures block.mean works
      r.ratio <- r.updated/W.r
      
      if (i <= n) {
        Q.star.updated <- force.sym(block.mean(Q.ratio, list(c(1:p, which(X[i, ] == 1)), p + q + 1)))
        r.star.updated <- c(mean(r.ratio[c(1:p, which(X[i, ] == 1))]), r.ratio[p + q + 1])
      } else {
        Q.star.updated <- force.sym(Q.ratio[c(p + i - n, p + q + 2), c(p + i - n, p + q + 2)])
        r.star.updated <- r.ratio[c(p + i - n, p + q + 2)]
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
    if (max.delta < abs.thresh && iteration > min.passes) {print("EP has converged; stopping EP"); break}
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
