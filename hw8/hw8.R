library(Matrix)
library(microbenchmark)
source("makeD2_sparse.R")
fmri = read.csv("fmri_z.csv")
colnames(fmri) = NULL
fmri = as.matrix(fmri)
c.fmri = c(fmri)

D = makeD2_sparse(nrow(fmri), ncol(fmri))
L = t(D) %*% D

lu.solve = function(A, b) {
  # sparseLU (see ?lu):
  # decomposition: A = P' L U Q
  # problem: A x = b
  # --------
  # P' L U Q x = b
  # P' L y = b
  # P P' L y = L y= P b ***
  # U Q x = y
  # U z = y ***
  # Q' Q x = x = Q' z ***
  # --------
  # So the results we need are
  # L y = P b
  # U z = y
  # x = Q' z
  lu.C = expand(lu(C))
  # Leverages a) sparse matrix functionality and b) L x = b solution
  # where L is lower-triangular
  y = Matrix::solve(lu.C$L, lu.C$P %*% c(fmri), system = "L")
  # Leverages a) sparse matrix functionality and b) L' x = b solution
  # where L' is upper-triangular
  z = Matrix::solve(lu.C$U, y, system = "Lt")

  x = t(lu.C$Q) %*% z
  return(x)
}

gaussseidel.solve = function(A, b, n.iter = 10) {
  # triu and tril take the upper and lower triangular matrices of their input
  # k = 0 by default keeps the diagonal, and k = 1 keeps only the diagonals
  # above the diagonal (k = 2 would throw out the diagonal above the main, etc.)
  U = triu(A, k = 1)
  L.star = tril(A)

  x = rep(0, length(b))
  for (i in 1:n.iter) {
    x.old = x
    # Leverages a) sparse matrix functionality and b) L x = b solution
    # where L is lower-triangular
    x = Matrix::solve(L.star, b - U %*% x, system = "L")
    if (sum(abs(x.old - x) ^ 2) < 1e-20) {
      #print(paste("gs:", i, "iterations"))
      break
    }
  }
  return(x)
}

jacobi.solve = function(A, b, n.iter = 100) {
  # Using the Diagonal constructor keeps everything sparse
  R = A - Diagonal(x = diag(A))
  # Only store the inverse of the values along the diagonal
  D.inv = 1 / diag(A)

  x = rep(0, length(b))
  for (i in 1:n.iter) {
    x.old = x
    x = D.inv * (b - R %*% x)
    if (sum(abs(x.old - x) ^ 2) < 1e-20) {
      #print(paste("jac:", i, "iterations"))
      break
    }
  }
  return(x)
}

lambda = 1
C = Diagonal(128 ** 2) + lambda * L
n.iter = 200
x = lu.solve(C, c(fmri))
gs.x = gaussseidel.solve(C, c(fmri), n.iter = n.iter)
jac.x = jacobi.solve(C, c(fmri), n.iter = n.iter)
# Plot
xMat = matrix(x, nrow = 128)
gsxMat = matrix(gs.x, nrow=128)
jacxMat = matrix(jac.x, nrow=128)

par(mfrow = c(2, 2))
cols = gray((0:15 / 15) ^ 0.67)
image(fmri, col = cols, main = "raw")
image(xMat, col = cols, main = "lu")
image(gsxMat, col = cols, main = "gs")
image(jacxMat, col = cols, main = "jac")

# HOBBITSES ARE TRICKSEY
# This is a completely unfair comparison as-is.
# the LU decomposition actually gets stored back into C (a property
# of the Matrix package). After lu(C) has been called,
# C@factors$LU caches the decomposition.
# I found this out after lu.solve was finishing in 8 ms and getting the right
# answer.
microbenchmark::microbenchmark(lu.solve(C, c.fmri),
                               gaussseidel.solve(C, c.fmri, n.iter = n.iter),
                               jacobi.solve(C, c.fmri, n.iter = n.iter),
                               times = 10L)

fair = function(fn, a, b, ...) {
  a@factors = list()
  fn(a, b, ...)
}
microbenchmark::microbenchmark(fair(lu.solve, C, c.fmri),
                               fair(gaussseidel.solve, C, c.fmri, n.iter = n.iter),
                               fair(jacobi.solve, C, c.fmri, n.iter = n.iter),
                               times = 10L)

