
# Auxiliary functions and variables for lasso linear regression

library(glmnet) # Fast ridge regression estimates

p.8 <- c(0.003246343272134, 0.051517477033972, 0.195077912673858, 0.315569823632818,
         0.274149576158423, 0.131076880695470, 0.027912418727972, 0.001449567805354)
s.8 <- c(1.365340806296348, 1.059523971016916, 0.830791313765644, 0.650732166639391,
         0.508135425366489, 0.396313345166341, 0.308904252267995, 0.238212616409306)

B.0 <- function(mu, sigma.2) {
  # Logistic-normal integral (0th order) approximation
  sum(p.8*pnorm(mu*s.8/sqrt(1 + sigma.2*s.8^2)))
}

B.1 <- function(mu, sigma.2) {
  # Logistic-normal integral (1st order) approximation
  sum(((p.8*sqrt(sigma.2)*s.8)/sqrt(1 + sigma.2*s.8^2))*dnorm(mu*s.8/sqrt(1 + sigma.2*s.8^2)))
}

B.2 <- function(mu, sigma.2) {
  # Logistic-normal integral (2nd order) approximation
  sum(p.8*(pnorm(mu*s.8/sqrt(1 + sigma.2*s.8^2)) - ((mu*sigma.2*s.8^3)/((1 + sigma.2*s.8^2)^1.5))*dnorm(mu*s.8/sqrt(1 + sigma.2*s.8^2))))
}

I.r <- function(y, m, V) {
  # Hybrid integrals for sites
  m <- as.numeric(m); v <- as.numeric(V)
  B.0 <- B.0(m, v)
  B.1 <- B.1(m, v)
  B.2 <- B.2(m, v)
  
  if (y == 1) {
    list(I.0 = B.0, 
         I.1 = m*B.0 + sqrt(v)*B.1, 
         I.2 = m^2*B.0 + 2*m*sqrt(v)*B.1 + v*B.2)
  } else {
    list(I.0 = 1 - B.0, 
         I.1 = m - (m*B.0 + sqrt(v)*B.1), 
         I.2 = m^2 + v - (m^2*B.0 + 2*m*sqrt(v)*B.1 + v*B.2))
  }
}

ep.approx <- function(X, y, sigma.2.beta, alpha, lambda.init,
                      max.passes, tol.factor, stop.factor, abs.thresh, 
                      rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian logistic regression
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  stop.ep <- F
  
  # Obtaining initial estimates
  init.mu <- as.vector(coef(glmnet(X[, -1], y, family = binomial(link = "logit"), alpha = 0, lambda = lambda.init)))
  init.Sigma <- diag(p)
  init.Sigma.inv <- force.sym(solve(init.Sigma))
  
  # Parameter initialisation
  Q.values <- array(dim = c(p, p, n + 1))
  r.values <- matrix(nrow = n + 1, ncol = p)
  
  for (i in 1:(n + 1)) {
    if (i <= n) {
      Q.values[, , i] <- (init.Sigma.inv - solve(sigma.2.beta*diag(p)))/n
      r.values[i, ] <- (init.Sigma.inv%*%init.mu)/n
    } else {
      Q.values[, , i] <- solve(sigma.2.beta*diag(p))
      r.values[i, ] <- rep(0, p)
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
      Q.cavity <- Q.sum - Q.values[, , i]
      Q.cavity.inv <- tryCatch(force.sym(solve(Q.cavity)), error = err)
      if (!is.matrix(Q.cavity.inv)) {stop.ep <- T; break}
      r.cavity <- r.sum - r.values[i, ]
      mu.cavity <- Q.cavity.inv%*%r.cavity
      Sigma.cavity <- Q.cavity.inv
      W <- X[i, ]
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
      I.r.res <- tryCatch(I.r(y[i], m, V), error = err)
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
      Q.updated <- Sigma.hybrid.inv - Q.cavity
      r.updated <- Sigma.hybrid.inv%*%mu.hybrid - r.cavity
      
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
