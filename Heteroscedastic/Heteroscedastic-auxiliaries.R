
# Auxiliary functions and variables for heteroscedastic linear regression

sigma.2.theta <- 10000

sim.settings <- list(c(n = 200, p.1 = 40, p.2 = 10),
                     c(n = 200, p.1 = 20, p.2 = 20),
                     c(n = 200, p.1 = 10, p.2 = 40))

bench.settings <- list(c(n = 40, p.1 = 2, p.2 = 2),
                       c(n = 725, p.1 = 6, p.2 = 2),
                       c(n = 125, p.1 = 5, p.2 = 5))

log.joint.likelihood <- function(theta, X.1, X.2, y, mu.theta, Sigma.theta) {
  # Log joint likelihood
  n <- nrow(X.1)
  p.1 <- ncol(X.1)
  p.2 <- ncol(X.2)
  
  beta.1 <- theta[1:p.1]
  beta.2 <- theta[(p.1 + 1):length(theta)]
  
  as.numeric(-0.5*(sum(2*X.2%*%beta.2) + t(y - X.1%*%beta.1)%*%((1/exp(2*X.2%*%beta.2))*(y - X.1%*%beta.1)) + t(theta - mu.theta)%*%solve(Sigma.theta)%*%(theta - mu.theta)))
}

laplace.approx <- function(X.1, X.2, y, mu.theta, Sigma.theta, lambda.init, maxit) {
  # Laplace approximation
  p.1 <- ncol(X.2)
  p.2 <- ncol(X.2)
  
  start <- if (p.1 > 2) c(as.vector(coef(glmnet(X.1[, -1], y, alpha = 0, lambda = lambda.init))), rep(0, p.2)) else
                        c(as.vector(coef(lm(y ~ X.1[, -1]))), rep(0, p.2))
  
  optim.res <- optim(start, fn = log.joint.likelihood, 
                     method = "BFGS", 
                     control = list(fnscale = -1, maxit = maxit),
                     hessian = T,
                     X.1 = X.1, X.2 = X.2, y = y, mu.theta = mu.theta, Sigma.theta = Sigma.theta)
  return(list(mu = optim.res$par, Sigma = solve(-optim.res$hessian)))
}

I.h.r <- function(y, mu, Sigma, eta, mult, lb.min, ub.max) {
  # Hybrid integrals for sites
  mu.1 <- mu[1]
  mu.2 <- mu[2]
  
  Q <- sym(solve(Sigma))
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

ep.approx <- function(X.1, X.2, y, mu.theta, Sigma.theta,
                      eta, alpha, Q.star.init, r.star.init, offset,
                      min.passes, max.passes, tol.factor, stop.factor,
                      abs.thresh, rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian heteroscedastic linear regression
  n <- nrow(X.1)
  p.1 <- ncol(X.1)
  p.2 <- ncol(X.2)
  stop.ep <- F
  pcount <- 0
  
  # Parameter initialisation
  Q.star.values <- array(dim = c(2, 2, n))
  r.star.values <- matrix(nrow = n, ncol = 2)
  
  Q.p <- sym(solve(Sigma.theta))
  r.p <- Q.p%*%mu.theta
  
  Q.dot <- sym(offset) + Q.p
  r.dot <- r.p
  
  for (i in 1:n) {
    A <- cbind(c(X.1[i, ], rep(0, p.2)), c(rep(0, p.1), X.2[i, ]))
    Q.star.values[, , i] <- sym(Q.star.init)
    r.star.values[i, ] <- r.star.init
    Q.dot <- Q.dot + sym(A%*%Q.star.values[, , i]%*%t(A))
    r.dot <- r.dot + A%*%r.star.values[i, ]
  }
  
  Sigma.dot <- sym(solve(Q.dot))
  mu.dot <- Sigma.dot%*%r.dot
  
  # Delta initialisation
  deltas <- matrix(-1, nrow = max.passes*n, ncol = 5)
  colnames(deltas) <- c("index", "iteration", "i", "delta", "skip")
  prev.max.delta <- Inf
  prev.med.delta <- Inf
  index <- 1
  
  # Main EP loop
  for (iteration in 1:max.passes) {
    
    if (verbose) print(paste0("---- Current iteration: ", iteration, " ----"))
    
    for (i in sample(1:n)) {
      
      A <- cbind(c(X.1[i, ], rep(0, p.2)), c(rep(0, p.1), X.2[i, ]))
      Sigma.dot.star <- sym(t(A)%*%Sigma.dot%*%A); mu.dot.star <- t(A)%*%mu.dot
      Q.dot.star <- sym(solve(Sigma.dot.star)); r.dot.star <- Q.dot.star%*%mu.dot.star
      Q.m.star <- Q.dot.star - eta*Q.star.values[, , i]; r.m.star <- r.dot.star - eta*r.star.values[i, ]
      Sigma.m.star <- sym(solve(Q.m.star)); mu.m.star <- Sigma.m.star%*%r.m.star
      
      I.h.r.res <- tryCatch(I.h.r(y[i], mu.m.star, Sigma.m.star, eta, mult = 5, lb.min = -5, ub.max = Inf), error = err)
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
