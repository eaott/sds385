library(Matrix)
source("laplacian_methods.R")
source("makeD2_sparse.R")

soft.threshold = function(x, lambda) {
  if ("numeric" %in% class(x)) {
    return(sign(x) * pmax(abs(x) - lambda, 0))
  }
  else if ("dgeMatrix" %in% class(x) && ncol(x) == 1) {
    return(sign(x) * pmax((abs(x) - lambda)[ , 1], 0))
  }
  stop(paste("Can only handle numeric vectors and dgCMatrix vectors", class(x)))
}

graphfusedlasso.solve = function(D, y, lambda = 1, n.iter = 200,
                                 tol.rel = 1e-5,
                                 tol.abs = 1e-5) {
  # Uses the distributed ADMM

  # Not sure about needing this:
  # if (!("dgCMatrix" %in% class(D))) {
  #   stop("Design matrix A should be a dgCMatrix")
  # }
  rho = 2 * lambda
  # minimize 1/2 * ||y - x||_2^2 + lambda ||Dx||_1
  x = rep(0, length(y))
  z = rep(0, length(y))
  u = rep(0, length(y))

  r = rep(0, nrow(D))
  s = rep(0, nrow(D))
  t = rep(0, nrow(D))

  # really need sparse D, and probably use an LU factorization instead
  # of the actual matrix.
  D.lu = Diagonal(ncol(D)) + t(D) %*% D

  residual.prim = c(x, r)
  residual.dual   = c(z, s)
  converged = FALSE
  iter = n.iter
  for (i in 1:n.iter) {
    z.prev = z
    s.prev = s

    x = 1 / (1 + rho) * (y + rho * (z - u))
    r = soft.threshold(s - t, lambda / rho)
    z = lu.solve(D.lu, x + u + t(D) %*% (r + t), useCache = TRUE)
    s = D %*% z
    u = u + x - z
    t = t + r - s

    ########################################
    # Convergence criteria
    ########################################
    residual.prim = c(x - z, (r - s)[ , 1])
    residual.dual = c(rho * (z - z.prev), rho * (s - s.prev)[ , 1])
    norm.prim = sum(residual.prim ^ 2)
    norm.dual = sum(residual.dual ^ 2)

    # Stopping conditions
    stop.prim = sqrt(length(y)) * tol.abs + tol.rel * max(c(x ^ 2, z ^ 2))
    stop.dual = sqrt(length(y)) * tol.abs + tol.rel * rho * sum(u ^ 2)
    # print(c(norm.prim, norm.dual, stop.prim, stop.dual))
    if (norm.prim <= stop.prim && norm.dual <= stop.dual) {
      converged = TRUE
      iter = i
      break
    }

    if (norm.prim >= 5 * norm.dual) {
      rho = 2 * rho
      u = u / 2
      t = t / 2
    }
    else if (norm.dual >= 5 * norm.prim) {
      rho = rho / 2
      u = u * 2
      t = t * 2
    }
  }
  if (!converged) {
    warning("Algorithm may not have converged. Try running for more iterations.")
  }
  return(list(x = x, z = z, r = r, s = s, n = iter))
}
