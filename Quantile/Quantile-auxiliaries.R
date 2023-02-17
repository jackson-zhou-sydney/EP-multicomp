
# EP code for Bayesian quantile regression

sourceCpp("Quantile/EP-approx.cpp")
sourceCpp("Quantile/MFVB-approx.cpp")

sigma.2.theta <- 10000
tau <- 0.75
sigma <- 0.2
n.binom <- 10

sim.settings <- list(c(n = 200, p = 40, dist = "normal"),
                     c(n = 200, p = 40, dist = "poisson"),
                     c(n = 200, p = 40, dist = "binomial"))

bench.settings <- list(c(n = 298, p = 2),
                       c(n = 235, p = 2),
                       c(n = 21, p = 4))

rho <- function(x, tau) {
  # Quantile loss function
  0.5*(abs(x) + (2*tau - 1)*x)
}

point.likelihood <- function(theta, x, y, tau) {
  # Likelihood evaluated at a point
  p <- length(x)
  
  beta <- theta[1:p]
  kappa <- theta[p + 1]
  
  as.numeric(tau*(1 - tau)*exp(-kappa - rho(y - t(x)%*%beta, tau)/exp(kappa)))
}

lppd <- function(X, y, tau, S) {
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
      subtotal <- subtotal + point.likelihood(S[s, ], x, y.point, tau)
    }
    
    total <- total + log(subtotal/n.s)
  }
  
  return(total)
}

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

I.h.r <- function(y, tau, mu, Sigma, eta, mult, lb.min, ub.max, length) {
  # Hybrid integrals for sites
  mu.1 <- mu[1]
  mu.2 <- mu[2]
  
  Q <- solve(Sigma)
  Q.11 <- Q[1, 1]
  Q.12 <- Q[1, 2]
  Q.22 <- Q[2, 2]
  
  abc <- function(x, lower) {
    c(Q.11,
      2*(Q.12*(x - mu.2) - Q.11*mu.1) + eta*(1 + ifelse(lower, -1, 1) - 2*tau)/exp(x),
      Q.11*mu.1^2 + 2*Q.12*mu.1*(mu.2 - x) + Q.22*(x - mu.2)^2 + eta*(2*x + (2*tau - 1 + ifelse(lower, 1, -1))*y/exp(x)))
  }
  
  lb <- max(mu.2 - mult*sqrt(Sigma[2, 2]), lb.min)
  ub <- min(mu.2 + mult*sqrt(Sigma[2, 2]), ub.max)
  
  x.values <- seq(from = lb, to = ub, length = length)
  y.matrix <- matrix(nrow = 6, ncol = length)
  
  for (i in 1:length(x.values)) {
    x <- x.values[i]
    abc.l <- abc(x, T)
    abc.u <- abc(x, F)
    y.matrix[1, i] <- TGI.minus.0(abc.l, y) + TGI.plus.0(abc.u, y)
    y.matrix[2, i] <- TGI.minus.1(abc.l, y) + TGI.plus.1(abc.u, y)
    y.matrix[3, i] <- x*y.matrix[1, i]
    y.matrix[4, i] <- TGI.minus.2(abc.l, y) + TGI.plus.2(abc.u, y)
    y.matrix[5, i] <- x*y.matrix[2, i]
    y.matrix[6, i] <- x^2*y.matrix[1, i]
  }
  
  ret.0 <- trapz(x.values, y.matrix[1, ])
  ret.11 <- trapz(x.values, y.matrix[2, ])
  ret.12 <- trapz(x.values, y.matrix[3, ])
  ret.211 <- trapz(x.values, y.matrix[4, ])
  ret.212 <- trapz(x.values, y.matrix[5, ])
  ret.222 <- trapz(x.values, y.matrix[6, ])
  
  list(I.h.0 = ret.0/sqrt(det(2*pi*Sigma)), 
       I.h.1 = c(ret.11, ret.12)/sqrt(det(2*pi*Sigma)), 
       I.h.2 = matrix(c(ret.211, ret.212, ret.212, ret.222)/sqrt(det(2*pi*Sigma)), nrow = 2))
}

ep.approx <- function(X, y, mu.theta, Sigma.theta, 
                      tau, eta, alpha, Q.star.init, r.star.init, offset,
                      min.passes, max.passes, tol.factor, stop.factor,
                      abs.thresh, rel.thresh, delta.limit, patience, verbose) {
  # Dampened power EP for Bayesian quantile regression
  n <- nrow(X)
  p <- ncol(X)
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
    A <- cbind(c(X[i, ], 0), c(rep(0, p), 1))
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
      
      A <- cbind(c(X[i, ], 0), c(rep(0, p), 1))
      Sigma.dot.star <- sym(t(A)%*%Sigma.dot%*%A); mu.dot.star <- t(A)%*%mu.dot
      Q.dot.star <- sym(solve(Sigma.dot.star)); r.dot.star <- Q.dot.star%*%mu.dot.star
      Q.m.star <- Q.dot.star - eta*Q.star.values[, , i]; r.m.star <- r.dot.star - eta*r.star.values[i, ]
      Sigma.m.star <- sym(solve(Q.m.star)); mu.m.star <- Sigma.m.star%*%r.m.star
      
      I.h.r.res <- tryCatch(I.h.r(y[i], tau, mu.m.star, Sigma.m.star, eta, mult = 5, lb.min = -10, ub.max = Inf, length = 50), error = err)
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
