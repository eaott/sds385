library(Matrix)
source("makeD2_sparse.R")

lu.solve = function(A, b, useCache = TRUE) {
  # Solves A x = b for x using a sparse LU decomposition.
  #
  # Args:
  #   A: Design matrix, should be a dgCMatrix for handling sparsity
  #   b: Response vector
  #   useCache: Determines whether to use any cached LU decomposition in A.
  #             In production, should be left alone, but it's useful for
  #             benchmarking to set to FALSE.
  #
  # Returns:
  #   The fitted value of x.
  if (!("dgCMatrix" %in% class(A))) {
    stop("Design matrix A should be a dgCMatrix")
  }

  # Hobbitses are tricksey
  # The lu() function will cache results in the associated Matrix object
  # and re-use them, which makes benchmarking a little unfair. This forces
  # the decomposition to happen every time.
  if (!useCache) {
    C@factors = list()
  }
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
  lu.A = expand(lu(A))
  # Leverages a) sparse matrix functionality and b) L x = b solution
  # where L is lower-triangular
  y = Matrix::solve(lu.A$L, lu.A$P %*% b, system = "L")
  # Leverages a) sparse matrix functionality and b) L' x = b solution
  # where L' is upper-triangular
  z = Matrix::solve(lu.A$U, y, system = "Lt")

  x = t(lu.A$Q) %*% z
  return(x)
}

gaussseidel.solve = function(A, b, n.iter = 200, tol = 1e-20) {
  # Solves A x = b for x using a the Gauss-Seidel iterative method
  #
  # Args:
  #   A: Design matrix, should be a dgCMatrix for handling sparsity
  #   b: Response vector
  #   n.iter: Max number of iterations to run.
  #   tol: Tolerance for sum of squared differences in guesses for x between
  #        iterations
  #
  # Returns:
  #   The fitted value of x.
  if (!("dgCMatrix" %in% class(A))) {
    stop("Design matrix A should be a dgCMatrix")
  }
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
    if (sum(abs(x.old - x) ^ 2) < tol) {
      #print(paste("gs:", i, "iterations"))
      break
    }
  }
  return(x)
}

jacobi.solve = function(A, b, n.iter = 200, tol = 1e-20) {
  # Solves A x = b for x using a the Jacobi iterative method
  #
  # Args:
  #   A: Design matrix, should be a dgCMatrix for handling sparsity
  #   b: Response vector
  #   n.iter: Max number of iterations to run.
  #   tol: Tolerance for sum of squared differences in guesses for x between
  #        iterations
  #
  # Returns:
  #   The fitted value of x.
  if (!("dgCMatrix" %in% class(A))) {
    stop("Design matrix A should be a dgCMatrix")
  }

  # Using the Diagonal constructor keeps everything sparse
  R = A - Diagonal(x = diag(A))
  # Only store the inverse of the values along the diagonal for faster
  # computation.
  D.inv = 1 / diag(A)

  x = rep(0, length(b))
  for (i in 1:n.iter) {
    x.old = x
    x = D.inv * (b - R %*% x)
    if (sum(abs(x.old - x) ^ 2) < tol) {
      #print(paste("jac:", i, "iterations"))
      break
    }
  }
  return(x)
}