########################################
# Section C
########################################
library(ggplot2)
library(microbenchmark)
library(reshape2)
library(Matrix)

source("matrix_inverse_functions.R")

# Run a check to see that a randomly generated matrix has the same solutions
# using the different methods.
N = 200
P = 100
X = matrix(rnorm(N * P), nrow = N)
W = diag(nrow = N)
y = rnorm(N)
# Ensure correctness
microbenchmark(inversion_method(X, y, W),
               smart_inversion_method(X, y, W),
               my_method(X, y, W),
               times = 1, control = list(order = "inorder"), check = my_check)

# This time, don't use the same matrices, just create them on the fly and
# compare the performance.
microbenchmark(compare(inversion_method, N, P),
               compare(smart_inversion_method, N, P),
               compare(my_method, N, P),
               times = 100)


# Now, generate some larger-scale benchmark data, getting the
# average time each algorithm took on an n x p matrix
# Note: this will take a while. Should be the case that
# min(ns) > max(ps) to ensure full rank.
results = compute_fairly_batch(c(3200),
                               c(10, 20, 40, 80, 160, 320, 640, 1280))
# Graph the data
names(results) = c("inverse", "smart inverse", "my method", "n", "p")
melted = melt(results, id.vars = c("n", "p"))
names(melted)[3] = "method"
ggplot(melted, aes(p, value, col = method)) +
  geom_line() + geom_point() +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10)) +
  scale_x_log10() +
  xlab("p") + ylab("average time (s)") + ggtitle("n=3200")

########################################
# Section D
########################################
# Showing how to create a sparse matrix
N = 2000
P = 1000
S = 0.05
X = matrix(rnorm(N * P), nrow = N)
X = X * matrix(rbinom(N * P, 1, S), nrow = N)
X = subset(melt(X), value != 0)
X = sparseMatrix(
  X$Var1,
  X$Var2,
  x = X$value)
X[1:5, 1:5]


results = compute_fairly_sparse_batch(c(1600),
                                      c(10, 20, 40, 80, 160, 320),
                                      c(0.05), iter = 1)

names(results) = c("inverse", "smart inverse", "my method",
                   "my sparse method", "n", "p", "s")

melted = melt(results, id.vars = c("n", "p", "s"))
names(melted)[4] = "method"
ggplot(melted, aes(p, value, col=method)) +
  geom_line() + geom_point() +
  scale_y_log10(breaks=c(0.001, 0.01, 0.1, 1, 10)) + scale_x_log10() +
  xlab("p") + ylab("average time (s)") + ggtitle("sparsity=0.05, n=1600")
