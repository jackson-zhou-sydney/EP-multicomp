
# Auxiliary functions and variables for heteroscedastic linear regression

library(glmnet) # Fast ridge regression estimates

log.joint.likelihood <- function(theta, X, Z, y, D) {
  # Log joint likelihood
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  beta <- theta[1:p]
  alpha <- theta[(p + 1):length(theta)]
  
  as.numeric(-0.5*(sum(Z%*%alpha) + t(y - X%*%beta)%*%((1/exp(Z%*%alpha))*(y - X%*%beta)) + t(theta)%*%D%*%theta))
}

ljl.grad <- function(theta, X, Z, y, D) {
  # Gradient of log joint likelihood
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  beta <- theta[1:p]
  alpha <- theta[(p + 1):length(theta)]
  
  ll.grad.beta <- -0.5*(-2*t(X)%*%((1/exp(Z%*%alpha))*(y - X%*%beta)))
  ll.grad.alpha <- -0.5*(rowSums(t(Z)) - t(Z*as.vector((y - X%*%beta)/exp(Z%*%alpha)))%*%(y - X%*%beta))
  
  as.vector(c(ll.grad.beta, ll.grad.alpha) - D%*%theta)
}

ljl.hess <- function(theta, X, Z, y, D) {
  # Hessian of log joint likelihood
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  
  beta <- theta[1:p]
  alpha <- theta[(p + 1):length(theta)]
  
  ll.hess.beta.beta <- -t(X*as.vector(1/exp(Z%*%alpha)))%*%X
  ll.hess.beta.alpha <- -t(X)%*%(Z*as.vector((y - X%*%beta)/exp(Z%*%alpha)))
  ll.hess.alpha.alpha <- -0.5*t(Z*as.vector(y - X%*%beta))%*%(Z*as.vector((y - X%*%beta)/exp(Z%*%alpha)))
  
  ll.hess <- matrix(nrow = p + q, ncol = p + q)
  ll.hess[1:p, 1:p] <- ll.hess.beta.beta
  ll.hess[1:p, (p + 1):(p + q)] <- ll.hess.beta.alpha
  ll.hess[(p + 1):(p + q), 1:p] <- t(ll.hess.beta.alpha)
  ll.hess[(p + 1):(p + q), (p + 1):(p + q)] <- ll.hess.alpha.alpha
  
  ll.hess - D
}

laplace.approx <- function(X, Z, y, D, start = rep(0, ncol(as.matrix(X)) + ncol(as.matrix(Z)))) {
  # Laplace approximation
  m <- optim(start, fn = log.joint.likelihood, gr = ljl.grad, 
             method = "BFGS", 
             control = list(fnscale = -1, maxit = 50000), 
             X = X, Z = Z, y = y, D = D)$par
  return(list(mu = m, Sigma = solve(-ljl.hess(m, X, Z, y, D))))
}

I.r <- function(y, m, V, eta, mult, maxEval, tol) {
  # Hybrid integrals for sites
  m.1 <- m[1]
  m.2 <- m[2]
  
  V.inv <- solve(V)
  V.11 <- V.inv[1, 1]
  V.12 <- V.inv[1, 2]
  V.22 <- V.inv[2, 2]
  
  abc <- function(x) {
    c(V.11 + eta/exp(x),
      2*(V.12*(x - m.2) - V.11*m.1) - eta*2*y/exp(x),
      V.11*m.1^2 + 2*V.12*m.1*(m.2 - x) + V.22*(x - m.2)^2 + eta*(y^2/exp(x) + x))
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

ep.approx <- function(X, Z, y, D, 
                      eta, alpha, lambda.init,
                      max.passes, tol.factor, stop.factor, abs.thresh, 
                      rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian heteroscedastic linear regression
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)
  stop.ep <- F
  
  # Obtaining initial estimates
  init.mu <- c(as.vector(coef(glmnet(X[, -1], y, alpha = 0, lambda = lambda.init))), rep(0, q))
  init.Sigma <- diag(p + q)
  init.Sigma.inv <- force.sym(solve(init.Sigma))
  
  # Parameter initialisation
  Q.values <- array(dim = c(p + q, p + q, n + 1))
  r.values <- matrix(nrow = n + 1, ncol = p + q)
  
  for (i in 1:(n + 1)) {
    if (i <= n) {
      Q.values[, , i] <- (init.Sigma.inv - D)/n
      r.values[i, ] <- (init.Sigma.inv%*%init.mu)/n
    } else {
      Q.values[, , i] <- D
      r.values[i, ] <- rep(0, p + q)
    }
  }
  
  Q.sum <- rowSums(Q.values, dims = 2)
  r.sum <- colSums(r.values)
  
  # Delta initialisation
  deltas <- matrix(-1, nrow = max.passes*n, ncol = 5)
  colnames(deltas) <- c("index", "iteration", "i", "delta", "skip")
  prev.max.delta <- Inf
  prev.med.delta <- Inf
  pcount <- 0
  index <- 1
  
  # Main EP loop
  for (iteration in 1:max.passes) {
    
    if (verbose) print(paste0("---- Current iteration: ", iteration, " ----"))
    
    for (i in sample(1:n)) {
      # Setting up
      Q.cavity <- Q.sum - eta*Q.values[, , i]
      Q.cavity.inv <- tryCatch(force.sym(solve(Q.cavity)), error = err)
      if (!is.matrix(Q.cavity.inv)) {stop.ep <- T; break}
      r.cavity <- r.sum - eta*r.values[i, ]
      mu.cavity <- Q.cavity.inv%*%r.cavity
      Sigma.cavity <- Q.cavity.inv
      W <- cbind(c(X[i, ], rep(0, q)), c(rep(0, p), Z[i, ]))
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
      I.r.res <- tryCatch(I.r(y[i], m, V, eta, mult = 5, maxEval = 0, tol = 0.0001), error = err)
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
      
      delta <- max(norm(r.updated - r.values[i, ], "2"), norm(Q.updated - Q.values[, , i], "F"))
      
      if (is.na(delta) || delta > tol.factor*prev.med.delta || delta > delta.limit) {
        deltas[index, ] <- c(index, iteration, i, delta, 1)
        index <- index + 1
        next
      } else {
        deltas[index, ] <- c(index, iteration, i, delta, 0)
        index <- index + 1
      }
      
      Q.new <- (1 - alpha)*Q.values[, , i] + alpha*Q.updated
      r.new <- (1 - alpha)*r.values[i, ] + alpha*r.updated
      Q.sum <- Q.sum - Q.values[, , i] + Q.new
      r.sum <- r.sum - r.values[i, ] + r.new
      Q.values[, , i] <- Q.new
      r.values[i, ] <- r.new
    }
    
    iteration.deltas <- deltas[deltas[, "iteration"] == iteration & deltas[, "skip"] == 0, ][, "delta"]
    max.delta <- max(iteration.deltas)
    med.delta <- median(iteration.deltas)
    if (verbose) {
      print(paste0("Maximum delta: ", round(max.delta, 2)))
      print(paste0("Median delta: ", round(med.delta, 2)))
    }
    
    if (stop.ep) {print("Too many numerical errors; stopping EP"); break}
    if (max.delta < abs.thresh) {print("EP has converged; stopping EP"); break}
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
