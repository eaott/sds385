wdbc = read.csv("wdbc.csv")

weight = function(x, b) {
  # x is NxP, b is Px1
  # x %*% b is Nx1
  w = (1/(1+exp(-as.matrix(x) %*% as.matrix(b))))
  w = pmax(w, 0.001)
  w = pmin(w, 0.999)
  return(w)
}

gradL = function(x, b, y, m) {
  # y is Nx1, m is 1x1, weight is Nx1
  # so this will multiply each element into the corresponding
  # row of x
  terms = as.numeric(y - m * weight(x, b)) * x
  # now, take the multiplied rows and sum them up creating
  # a 1xP vector that is just a c(...) so Px1
  return(-colSums(terms))
}

logL = function(x, y, b, m) {
  # probably do sum of log pmfs 
  # with dbinom(y, m, weight(x, b))
  return(sum(dbinom(y, m, weight(x, b), log = TRUE)))
}

X = wdbc[ , 3:12]
# Add a column of 1s.
X = unname(as.matrix(cbind(1, X)))
Y = wdbc[ , 2]
# Y: replace M with 1 and B with 0.
Y = as.matrix((Y == "M") * 1)

N = dim(X)[1]
P = dim(X)[2]
beta = as.matrix(rnorm(P))
iter = 100
logLike = rep(0, iter)
for (i in 1:iter) {
  # check for +/-
  beta = beta + 0.0001 * gradL(X, beta, Y, 1)
  logLike[i] = logL(X, Y, beta, 1)
}

plot(logLike)

