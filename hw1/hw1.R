library(microbenchmark)

inversion_method = function(X, y, W) {
  A = solve(t(X) %*% W %*% X)
  return(A %*% t(X) %*% W %*% y)
}

smart_inversion_method = function(X, y, W) {
  # Uses just the diagonal of W to not multiply all the zeroes.
  w = diag(W)
  A = solve(t(X) %*% (w * X))
  return(A %*% t(X) %*% (w * y))
}

my_method = function(X, y, W) {
  # TODO: update
  # from http://www.seas.ucla.edu/~vandenbe/103/lectures/chol.pdf
  w = diag(W)
  A = t(X) %*% (w * X)
  R = chol(A)
  L = t(R)
  z = forwardsolve(L, t(X) %*% (w * y))
  return(backsolve(R, z))
}


my_check <- function(values) {
  # In cases where the values should all be near-identical
  bhat = values[[1]]
  all(sapply(values[-1], function(x) {sum((x-bhat)^2) < 1e-9}))
}

compare = function(method, n, p) {
  X = matrix(rnorm(N * P), nrow=N)
  W = diag(nrow=N)
  y = rnorm(N)
  return(method(X, y, W))
}
N = 2000
P = 1000
X = matrix(rnorm(N * P), nrow=N)
W = diag(nrow=N)
y = rnorm(N)
# Ensure correctness
microbenchmark(inversion_method(X, y, W), smart_inversion_method(X, y, W), my_method(X, y, W), times=1, control=list(order="inorder"), check=my_check)

# Better timing
microbenchmark(compare(inversion_method, N, P), compare(smart_inversion_method, N, P), compare(my_method, N, P), times=100)

