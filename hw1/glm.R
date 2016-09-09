
weight = function(x, b) {
  # x is NxP, b is Px1
  # x %*% b is Nx1
  w = (1/(1+exp(-as.matrix(x) %*% as.matrix(b))))
  w = pmax(w, 1e-6)
  w = pmin(w, 1 - 1e-6)
  return(w)
}


gradL = function(x, b, y, m) {
  # y is Nx1, m is 1x1, weight is Nx1
  # so this will multiply each element into the corresponding
  # row of x
  terms = as.numeric(y - m * weight(x, b)) * x
  # now, take the multiplied rows and sum them up creating
  # a 1xP vector that is just a c(...) so Px1
  gL = colSums(terms)
  return(-gL)
}

gradLNewton = function(x, b, y, m) {
  w = c(weight(x, b))
  hessian = t(x) %*% (m * w * (1-w) * x)
  gL = as.matrix(gradL(x, b, y, m))
  # FIXME want solve(hessian, gL)
  return(-solve(hessian) %*% gL)
}


logL = function(x, b, y, m) {
  # probably do sum of log pmfs
  # with dbinom(y, m, weight(x, b))
  return(sum(dbinom(y, m, weight(x, b), log = TRUE)))
}

wdbc = read.csv("wdbc.csv")

X = scale(wdbc[ , 3:12])
#X = wdbc[ , 3:12]

# Add a column of 1s.
X = unname(as.matrix(cbind(1, X)))
Y = wdbc[ , 2]
# Y: replace M with 1 and B with 0.
Y = as.matrix((Y == "M") * 1)

N = dim(X)[1]
P = dim(X)[2]
newton = FALSE
iter = 20000
allYourBeta = c()
for (j in 1:1) {
  beta = rep(0, P)
  logLike = rep(0, iter)
  for (i in 1:iter) {
    if (newton) {
      beta = beta + gradLNewton(X, beta, Y, 1)
    } else {
      gL = gradL(X, beta, Y, 1)
      # Scaled version
      beta = beta - 0.001 * gL
      # Unscaled
      # beta = beta - 0.001 * gL / sqrt(sum(gL ^ 2))
    }
    logLike[i] = logL(X, beta, Y, 1)
  }
  plot(logLike, type="l")
  allYourBeta = cbind(allYourBeta, as.numeric(beta))
}

# Scaled X, gradient
# 20000 iters
# LogL = -74.10653
#  [1] -0.2504395 -0.5149440  1.6259505 -0.9047349  5.4986217
# [6]  1.0729947 -0.5458339  0.9790104  2.3175757  0.4692481
# [11] -0.2735361

# Scaled X, newton
# 200 iters
# LogL = -73.06523
# [1]  0.46964720 -7.22158870  1.64995454 -1.73632394 14.00606994
# [6]  1.07358670 -0.07657184  0.67150243  2.58049248  0.44471461
# [11] -0.48075207
