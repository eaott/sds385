library(Rcpp)
library(RcppEigen)
library(Matrix)
########################################
# Using the BIG data
########################################

X = t(readRDS('url_x.rds'))
y = readRDS('url_y.rds')

Rcpp::sourceCpp("full_scale.cpp")
t = proc.time()
ret = test(X, y)
proc.time() - t


#ret = adagrad_sparse_sgd(X, y, rep(0, ncol(X)), rep(0, ncol(X)), 0.1)

