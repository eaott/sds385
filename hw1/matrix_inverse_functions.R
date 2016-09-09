library(ggplot2)
library(microbenchmark)
library(reshape2)
library(Matrix)


inversion_method = function(X, y, W) {
  # Uses the raw inverse, which is more unstable.
  #
  # Args:
  #   X: Design matrix, N x P
  #   y: Response vector, N x 1
  #   W: Weight matrix, N x N
  #
  # Returns:
  #   The solution to the system of equations
  A = solve(t(X) %*% W %*% X)
  return(A %*% t(X) %*% W %*% y)
}

smart_inversion_method = function(X, y, W) {
  # Uses just the diagonal of W to not multiply all the zeroes.
  #
  # Args:
  #   X: Design matrix, N x P
  #   y: Response vector, N x 1
  #   W: Weight matrix, N x N
  #
  # Returns:
  #   The solution to the system of equations
  w = diag(W)
  A = solve(t(X) %*% (w * X))
  return(A %*% t(X) %*% (w * y))
}

my_method = function(X, y, W) {
  # A smarter method that leverages the Cholesky matrix decomposition
  # instead of the traditional inverse.
  #
  # Args:
  #   X: Design matrix, N x P
  #   y: Response vector, N x 1
  #   W: Weight matrix, N x N
  #
  # Returns:
  #   The solution to the system of equations.
  # Notes:
  #   from http://www.seas.ucla.edu/~vandenbe/103/lectures/chol.pdf
  w = diag(W)
  A = t(X) %*% (w * X)
  R = chol(A)
  L = t(R)
  z = forwardsolve(L, t(X) %*% (w * y))
  return(backsolve(R, z))
}


my_check <- function(values) {
  # Function to check whether results in microbenchmark are the same.
  #
  # Args:
  #   values: list of results from microbenchmark. This function
  #           should not be called directly by the user
  #
  # Returns:
  #   TRUE if the sum of squared differences from the first item are within
  #   a tolerance of 1e-9.
  bhat = values[[1]]
  all(sapply(values[-1], function(x) {sum((x - bhat) ^ 2) < 1e-9}))
}

get_time = function() {
  # Simple function for determining execution time in seconds.
  #
  # Returns:
  #   Current elapsed time.
  return(unname(proc.time()["elapsed"]))
}

compare = function(method, n, p) {
  # Wrapper function that creates an n x p design matrix and calls the
  # given function and returns its result.
  #
  # Args:
  #   method: function that takes (x, y, w) arguments, where x is n x p,
  #           y is n x 1, and w is n x n
  #   n: number of rows
  #   p: number of columns
  #
  # Returns:
  #   The solution to the system of equations generated.
  X = matrix(rnorm(n * p), nrow = n)
  W = diag(nrow = n)
  y = rnorm(n)
  return(method(X, y, W))
}

compare_fairly = function(n, p, iter=100) {
  # Replacement for microbenchmark that provides for re-using the same matrix.
  #
  # Args:
  #   n: number of rows
  #   p: number of columns
  #   iter: number of iterations (examples to generate)
  #
  # Returns:
  #   data frame indicating how fast each algorithm was on average: $itime is
  #   the basic inverse, $stime is the "smart" inverse that uses the diagonal
  #   matrix in an efficient way, and $mtime is my new method (Cholesky)
  itime = 0
  stime = 0
  mtime = 0
  for (i in 1:iter) {
    X = matrix(rnorm(n * p), nrow = n)
    W = diag(nrow = n)
    y = rnorm(n)

    t = get_time()
    inversion_method(X, y, W)
    itime = itime + (get_time() - t)

    t = get_time()
    smart_inversion_method(X, y, W)
    stime = stime + (get_time() - t)

    t = get_time()
    my_method(X, y, W)
    mtime = mtime + (get_time() - t)
  }
  return(data.frame(itime = c(itime / iter),
              stime = c(stime / iter),
              mtime = c(mtime / iter)))
}

compute_fairly_batch = function(ns, ps, iter=100) {
  # Essentially does the same as compare_fairly, but for multiple sizes
  # of matrices
  #
  # Args:
  #   ns: values of n to try
  #   ps: values of p to try
  #   iter: number of iterations to use in each attempt
  #
  # Returns:
  #   Results from compare_fairly in a big data frame indicating the size
  #   of the matrix
  results = data.frame()
  for (n in ns) {
    for (p in ps) {
      results = rbind(results, cbind(compare_fairly(n, p, iter = iter),
                                     list(n = n, p = p)))
      print(paste(n, p))
    }
  }
  return(results)
}


########################################
# Section D
########################################
my_method_sparse = function(X, y, W) {
  # Uses a sparse matrix representation with the Cholesky decomposition.
  #
  # Args:
  #   X: Design matrix, N x P
  #   y: Response vector, N x 1
  #   W: Weight matrix, N x N
  #
  # Returns:
  #   Solution to the system of equations
  myDims = dim(X)
  X = subset(melt(X), value != 0)
  X = Matrix::sparseMatrix(
    X$Var1,
    X$Var2,
    x = X$value,
    dims = myDims)
  # from http://www.seas.ucla.edu/~vandenbe/103/lectures/chol.pdf
  # Need t(X) %*% W %*% X, so break up into
  # [t(X) %*% t(W^(1/2))] %*% [W^(1/2) %*% X]
  w = diag(W) ^ 0.5
  # does t(A) %*% A faster
  A = crossprod(w * X)
  R = chol(A)
  L = t(R)
  z = forwardsolve(L, t(X) %*% (w * y))
  return(backsolve(R, z))
}

compare_fairly_sparse = function(n, p, s, iter=100) {
  # Same as compare_fairly, but also includes the sparse matrix
  #
  # Args:
  #   n: num rows
  #   p: num cols
  #   s: sparsity (1 = fully dense, 0 = all zeroes)
  #   iter: number of examples to generate
  #
  # Returns:
  #   Returns average time for all methods. Same as compare_fairly, but includes
  #   $sparsetime for the new sparse method
  itime = 0
  stime = 0
  mtime = 0
  sparsetime = 0
  for (i in 1:iter) {
    X = matrix(rnorm(n * p), nrow = n)
    W = diag(nrow = n)
    y = rnorm(n)
    X = X * matrix(rbinom(n * p, 1, s), nrow = n)

    t = get_time()
    inversion_method(X, y, W)
    itime = itime + (get_time() - t)

    t = get_time()
    smart_inversion_method(X, y, W)
    stime = stime + (get_time() - t)

    t = get_time()
    my_method(X, y, W)
    mtime = mtime + (get_time() - t)

    t = get_time()
    my_method_sparse(X, y, W)
    sparsetime = sparsetime + (get_time() - t)
  }
  return(data.frame(itime = c(itime / iter),
                    stime = c(stime / iter),
                    mtime = c(mtime / iter),
                    sparsetime = c(sparsetime / iter)))
}

compute_fairly_sparse_batch = function(ns, ps, ss, iter=100) {
  # Same as compute_fairly_batch, but also includes sparsity.
  #
  # Args:
  #   ns: sets of n
  #   ps: sets of p
  #   ss: sets of s (sparsity)
  #   iter: number of iterations
  #
  # Returns:
  #   Returns results from compute_fairly_sparse, but along with n, p, and s.
  results = data.frame()
  for (n in ns) {
    for (p in ps) {
      for (s in ss) {
        results = rbind(results,
                        cbind(compare_fairly_sparse(n, p, s),
                              list(n = n, p = p, s = s)))
        print(paste(n, p, s))
      }
    }
  }
  return(results)
}
