wdbc = read.csv("wdbc.csv")

weight = function(x, b) {
  # x is NxP, b is Px1
  # x %*% b is Nx1
  w = (1/(1+exp(-as.matrix(x) %*% as.matrix(b))))
  w = pmax(w, 1e-6)
  w = pmin(w, 1 - 1e-6)
  return(w)
}

gradL = function(x, b, y, m) {
  # FIXME something is still wrong here
  
  # y is Nx1, m is 1x1, weight is Nx1
  # so this will multiply each element into the corresponding
  # row of x
  terms = as.numeric(y - m * weight(x, b)) * x
  # now, take the multiplied rows and sum them up creating
  # a 1xP vector that is just a c(...) so Px1
  return(-colSums(terms))
}

logL = function(x, y, b, m) {
  # FIXME this is also broken, giving -Inf.
  
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
allYourBeta = c()
for (j in 1:10) {
  beta = as.matrix(1/colMeans(X) * sample(c(-1, 1), P, replace=TRUE))
  iter = 10000
  logLike = rep(0, iter)
  for (i in 1:iter) {
    # check for +/-?
    beta = beta - 1e-7 * gradL(X, beta, Y, 1)
    logLike[i] = logL(X, Y, beta, 1)
  }
  plot(logLike, type="l")
  allYourBeta = cbind(allYourBeta, as.numeric(beta))
}


