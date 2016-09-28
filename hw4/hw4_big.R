library(Rcpp)
library(RcppEigen)
library(Matrix)
########################################
# Using the BIG data
########################################

X = readRDS('url_x.rds')
y = readRDS('url_y.rds')
tX = t(X)

Rcpp::sourceCpp("full_scale.cpp")
t = proc.time()
test(tX)
proc.time() - t
t = proc.time()
testFn(tX)
proc.time() - t

ret = adagrad_sparse_sgd(X, y, rep(0, ncol(X)), rep(0, ncol(X)), 0.1)

