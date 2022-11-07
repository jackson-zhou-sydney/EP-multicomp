
# EP code for Bayesian quantile regression

sigma <- 0.2
n.binom <- 10

sim.settings <- list(c(n = 200, p = 40, dist = "normal"),
                     c(n = 200, p = 40, dist = "poisson"),
                     c(n = 200, p = 40, dist = "binomial"))

bench.settings <- list(c(n = 442, p = 11),
                       c(n = 97, p = 9),
                       c(n = 120, p = 201))

expec.lnig <- function(A, B, C, D, E, fun, radius = 100, lb) {
  # Expectation of product of reparameterised log-normal and inverse gamma densities
  p <- function(x) -A/exp(2*x) + B/exp(x) - C*x - (x - D)^2/(2*E)
  p.max <- optimise(p, c(-radius, radius), maximum = TRUE)$objective
  q <- function(x) exp(-A/exp(2*x) + B/exp(x) - C*x - (x - D)^2/(2*E) - p.max)
  
  if (fun == "x") {
    r <- function(x) exp(-A/exp(2*x) + B/exp(x) - C*x - (x - D)^2/(2*E) - p.max)*x
  } else if (fun == "x^2") {
    r <- function(x) exp(-A/exp(2*x) + B/exp(x) - C*x - (x - D)^2/(2*E) - p.max)*x^2
  } else if (fun == "1/exp(x)") {
    r <- function(x) exp(-A/exp(2*x) + B/exp(x) - C*x - (x - D)^2/(2*E) - p.max - x)
  } else if (fun == "1/exp(2*x)") {
    r <- function(x) exp(-A/exp(2*x) + B/exp(x) - C*x - (x - D)^2/(2*E) - p.max - 2*x)
  } else {
    stop("fun must be one of: x, x^2, 1/exp(x) or 1/exp(2*x)")
  }
  
  integrate(r, lb, Inf, abs.tol = 0)$value/integrate(q, lb, Inf, abs.tol = 0)$value
}

mfvb.approx <- function(X, y, mu.beta, Sigma.beta, mu.kappa, sigma.2.kappa,
                        tau, maxit, tol, verbose) {
  # MFVB for Bayesian quantile regression
  n <- nrow(X)
  p <- ncol(X)
  
  mu.beta.p <- mu.beta
  Sigma.beta.p <- Sigma.beta
  Sigma.beta.p.inv <- solve(Sigma.beta.p, tol = 1.0E-99)
  X.colSums <- colSums(X)
  
  # Initialisations
  mu.beta <- rep(0, p)
  Sigma.beta <- diag(p)
  expec.ie.kappa <- 1
  expec.ie.2.kappa <- 1
  expec.a <- rep(1, n)
  
  expec.y.X.beta.2 <- as.vector((y - X%*%mu.beta)^2 + diag(X%*%Sigma.beta%*%t(X)))
  A <- diag(expec.a)
  
  # Main MFVB loop
  for (iteration in 1:maxit) {
    # Storing old values
    mu.beta.old <- mu.beta
    
    # Update q(beta)
    Q.inv <- expec.ie.2.kappa*tau*(1 - tau)*t(X)%*%A%*%X + Sigma.beta.p.inv
    Q <- solve(Q.inv, tol = 1.0E-99)
    mu.beta <- as.vector(Q%*%(expec.ie.2.kappa*tau*(1 - tau)*t(X)%*%A%*%y - 
                              expec.ie.kappa*(0.5 - tau)*X.colSums + 
                              Sigma.beta.p.inv%*%mu.beta.p))
    Sigma.beta <- Q
    expec.y.X.beta.2 <- as.vector((y - X%*%mu.beta)^2 + diag(X%*%Sigma.beta%*%t(X)))
    
    # Update q(kappa)
    expec.ie.kappa <- expec.lnig(0.5*tau*(1 - tau)*as.numeric(t(expec.a)%*%expec.y.X.beta.2),
                                 (0.5 - tau)*sum(y - X%*%mu.beta), n, mu.kappa, sigma.2.kappa,
                                 fun = "1/exp(x)", lb = -5)
    expec.ie.2.kappa <- expec.lnig(0.5*tau*(1 - tau)*as.numeric(t(expec.a)%*%expec.y.X.beta.2),
                                   (0.5 - tau)*sum(y - X%*%mu.beta), n, mu.kappa, sigma.2.kappa,
                                   fun = "1/exp(2*x)", lb = -5)
    
    # Update q(a)
    expec.a <- 0.5/(tau*(1 - tau))*(expec.ie.2.kappa*expec.y.X.beta.2)^-0.5
    A <- diag(expec.a)
    
    # Checking for convergence
    delta <- norm(mu.beta - mu.beta.old, "2")
    if (verbose) print(paste0("Iteration ", iteration, ", current delta is: ", delta))
    if (delta < tol) break
  }
  
  # Returning parameters
  mu.kappa.q <- expec.lnig(0.5*tau*(1 - tau)*as.numeric(t(expec.a)%*%expec.y.X.beta.2),
                           (0.5 - tau)*sum(y - X%*%mu.beta), n, mu.kappa, sigma.2.kappa,
                           fun = "x", lb = -5)
  sigma.2.kappa.q <- expec.lnig(0.5*tau*(1 - tau)*as.numeric(t(expec.a)%*%expec.y.X.beta.2),
                                (0.5 - tau)*sum(y - X%*%mu.beta), n, mu.kappa, sigma.2.kappa,
                                fun = "x^2", lb = -5) - mu.kappa.q^2
  
  mu.theta <- numeric(p + 1)
  mu.theta[1:p] <- mu.beta
  mu.theta[p + 1] <- mu.kappa.q
  
  Sigma.theta <- matrix(0, nrow = p + 1, ncol = p + 1)
  Sigma.theta[1:p, 1:p] <- Sigma.beta
  Sigma.theta[p + 1, p + 1] <- sigma.2.kappa.q
  
  return(list(mu = mu.theta, Sigma = Sigma.theta))
}

I.r <- function(y, tau, m, V, eta, mult, lb.min, ub.max) {
  # Hybrid integrals for sites
  m.1 <- m[1]
  m.2 <- m[2]
  
  V.inv <- solve(V)
  V.11 <- V.inv[1, 1]
  V.12 <- V.inv[1, 2]
  V.22 <- V.inv[2, 2]
  
  abc <- function(x, lower) {
    c(V.11,
      2*(V.12*(x - m.2) - V.11*m.1) + eta*(1 + ifelse(lower, -1, 1) - 2*tau)/exp(x),
      V.11*m.1^2 + 2*V.12*m.1*(m.2 - x) + V.22*(x - m.2)^2 + eta*(2*x + (2*tau - 1 + ifelse(lower, 1, -1))*y/exp(x)))
  }
  
  lb <- max(m.2 - mult*sqrt(V[2, 2]), lb.min)
  ub <- min(m.2 + mult*sqrt(V[2, 2]), ub.max)
  
  ret.0 <- integrate(Vectorize(function(x) TGI.minus.0(abc(x, T), y) + TGI.plus.0(abc(x, F), y)), lb, ub, abs.tol = 0)$value
  ret.11 <- integrate(Vectorize(function(x) TGI.minus.1(abc(x, T), y) + TGI.plus.1(abc(x, F), y)), lb, ub, abs.tol = 0)$value
  ret.12 <- integrate(Vectorize(function(x) x*(TGI.minus.0(abc(x, T), y) + TGI.plus.0(abc(x, F), y))), lb, ub, abs.tol = 0)$value
  ret.211 <- integrate(Vectorize(function(x) TGI.minus.2(abc(x, T), y) + TGI.plus.2(abc(x, F), y)), lb, ub, abs.tol = 0)$value
  ret.212 <- integrate(Vectorize(function(x) x*(TGI.minus.1(abc(x, T), y) + TGI.plus.1(abc(x, F), y))), lb, ub, abs.tol = 0)$value
  ret.222 <- integrate(Vectorize(function(x) x^2*(TGI.minus.0(abc(x, T), y) + TGI.plus.0(abc(x, F), y))), lb, ub, abs.tol = 0)$value
  
  list(I.0 = ret.0/sqrt(det(2*pi*V)), 
       I.1 = c(ret.11, ret.12)/sqrt(det(2*pi*V)), 
       I.2 = matrix(c(ret.211, ret.212, ret.212, ret.222)/sqrt(det(2*pi*V)), nrow = 2))
}

ep.approx <- function(X, y, mu.theta, Sigma.theta, 
                      tau, eta, alpha, Q.star, r.star, prec,
                      min.passes, max.passes, tol.factor, stop.factor, 
                      abs.thresh, rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian quantile regression
  n <- nrow(X)
  p <- ncol(X)
  stop.ep <- F
  
  # Parameter initialisation
  Q.star.values <- array(dim = c(2, 2, n))
  r.star.values <- matrix(nrow = n, ncol = 2)
  
  Q.p <- force.sym(solve(Sigma.theta))
  r.p <- Q.p%*%mu.theta
  
  Q.sum <- prec*diag(p + 1) + Q.p
  r.sum <- r.p
  
  for (i in 1:n) {
    W <- cbind(c(X[i, ], 0), c(rep(0, p), 1))
    Q.star.values[, , i] <- Q.star
    r.star.values[i, ] <- r.star
    Q.sum <- Q.sum + force.sym(W%*%Q.star.values[, , i]%*%t(W))
    r.sum <- r.sum + W%*%r.star.values[i, ]
  }
  
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
      W <- cbind(c(X[i, ], 0), c(rep(0, p), 1))
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
      I.r.res <- tryCatch(I.r(y[i], tau, m, V, eta, mult = 5, lb.min = -10, ub.max = Inf), error = err)
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
      Q.ratio <- Q.updated/(W.r%*%t(W.r)); Q.ratio[!is.finite(Q.ratio)] <- NA
      r.ratio <- r.updated/W.r; r.ratio[!is.finite(r.ratio)] <- NA
      
      Q.star.updated <- force.sym(block.mean(Q.ratio, list(1:p, p + 1), na.rm = T))
      r.star.updated <- c(mean(r.ratio[1:p], na.rm = T), r.ratio[p + 1])
      
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
