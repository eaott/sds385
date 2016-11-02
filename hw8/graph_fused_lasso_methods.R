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

graphfusedlasso.solve = function(D, y, lambda = 1, n.iter = 200, tol = 1e-20) {
  # Uses the distributed ADMM

  # Not sure about needing this:
  # if (!("dgCMatrix" %in% class(D))) {
  #   stop("Design matrix A should be a dgCMatrix")
  # }
  rho = 1
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
  for (i in 1:n.iter) {
    x = 1 / (1 + rho) * (y + rho * (z - u))
    r = soft.threshold(s - t, lambda / rho)
    z = lu.solve(D.lu, (x + u + t(D) %*% (r + t)), useCache = TRUE)
    s = D %*% z
    u = u + x - z
    t = t + r - s

    ########################################
    # FIXME: convergence criteria
    ########################################

  }
  return(x)
}
