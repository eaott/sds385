findD = function(x, lambda, TOL = 1e-20, TOL.d = NULL) {
  if ("dgeMatrix" %in% class(x)) {
    x = x[ , 1]
  }
  l2.x = sqrt(sum(x ^ 2))
  if (l2.x == 0) {
    stop("Bad x vector")
  }

  # in this case, S(x, 0) = x
  if (sum(abs(x / l2.x)) <= lambda) {
    return(x / l2.x)
  }

  # if that's not the case, we binary search
  dmax = max(abs(x))
  dmin = 0
  if (is.null(TOL.d)) {
    TOL.d = 1e-6 * dmax
  }
  # scoping
  u = rep(0, length(x))
  while (dmax - dmin > TOL.d) {
    # find midpoint of range
    d = (dmax + dmin) / 2

    ########################################
    # FIXME: this is the offending line in S4
    ########################################

    tempU = sign(x) * pmax(abs(x) - d, 0)
    u = tempU / sqrt(sum(tempU ^ 2))
    l1.u = sum(abs(u))
    if (abs(l1.u - lambda) < TOL) {
      break
    }
    else if (l1.u < lambda) {
      # d needs to relax -- get closer to zero
      dmax = d
    }
    else if (l1.u > lambda) {
      # d needs to get stronger -- closer to full
      dmin = d
    }
  }

  return(u)
}

pmd.1 = function(X, lambda.u, lambda.v, n.iter = 10, TOL = 1e-10) {
  N = nrow(X)
  P = ncol(X)
  v = rnorm(P)
  v = v / sqrt(sum(v ^ 2))
  u = rep(0, N)
  i = 1
  for (i in 1:n.iter) {
    old.u = u
    old.v = v
    Xv = X %*% v

    ########################################
    # FIXME: squirrely business with d1
    #
    # where d1 = 0 if \norm{u}_1 \leq c1
    # otherwise, d1 is a positive constant such that \norm{u}_1 = c1
    #
    # find by binary search:
    # know d1 >= 0
    #
    # if d1 > max(abs(Xv)) then abs(Xv) - d1 < 0 so tempU = 0
    ########################################
    u = findD(Xv, lambda.u)

    Xu = t(X) %*% u
    v = findD(Xu, lambda.v)

    if (sum(abs(old.u - u)) < TOL && sum(abs(old.v - v)) < TOL) {
      break
    }
  }
  if (i == n.iter) {
    warning("may not have converged")
  }
  d = t(u) %*% X %*% v
  return(list(u = u,
              v = v,
              d = d[1, 1]))
}

pmd = function(X, lambda.u, lambda.v, rank, n.iter = 10, TOL = 1e-10) {
  u = c()
  v = c()
  d = c()
  for (i in 1:rank) {
    res = pmd.1(X, lambda.u, lambda.v, n.iter = n.iter, TOL = TOL)
    u = cbind(u, res$u)
    v = cbind(v, res$v)
    d = c(d, res$d)
    X = X - res$u %*% t(res$v) * res$d
  }
  return(list(u = u,
              v = v,
              d = d,
              resid = X))
}