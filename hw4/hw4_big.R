library(Rcpp)
library(RcppEigen)
library(Matrix)
########################################
# Using the BIG data
########################################

X = t(readRDS('url_x.rds'))
X = rbind(1, X)
y = readRDS('url_y.rds')

# Keep 20% for testing.
training = 1 - (runif(ncol(X)) < 0.2)

Rcpp::sourceCpp("full_scale.cpp")
nIter = 1
t = proc.time()
ret = sgd(X, y, training, nIter, 0.1)
proc.time() - t

ll = rep(0, length(ret$logL))
ll[1] = ret$logL[1]
for (i in 2:length(ll)) {
  #ll[i] = ((i - 1) * ll[i - 1] +  ret$logL[i]) / i
  ll[i] = 0.99 * ll[i - 1] + 0.01 * ret$logL[i]
}
plot(-ll * length(ll) / nIter, type="l")
tail(-ll * length(ll) / nIter)

########################################
# FIXME: Can add an accuracy calculation at the end
# computing \hat{y} < 0.5 and using y.
########################################








########################################
# Use the old data to show that we get similar
# results to glm.
########################################
wdbc = read.csv("wdbc.csv")

# Switch between using a scaled and non-scaled design matrix
wX = scale(wdbc[ , 3:12])
#X = wdbc[ , 3:12]
# Add a column of 1s.
wX = unname(as.matrix(cbind(1, wX)))

wY = wdbc[ , 2]
# Y: replace M with 1 and B with 0.
wY = (wY == "M") * 1
wX = t(wX)
wdcgX = Matrix(wX, sparse=TRUE)

nIter = 10000
wRet = sgd(wdcgX, wY, rep(1, ncol(wX)), nIter, 0.1)

ll = rep(0, length(wRet$logL))
ll[1] = wRet$logL[1]
for (i in 2:length(ll)) {
  #ll[i] = ((i - 1) * ll[i - 1] +  wRet$logL[i]) / i
  ll[i] = 0.99 * ll[i - 1] + 0.01 * wRet$logL[i]
}
plot(-ll * length(ll) / nIter, type="l")
tail(-ll * length(ll) / nIter)
logL(t(wX), wRet$beta, wY, 1)
logL(t(wX), coef(glm(wY ~ t(wX) - 1, family="binomial")), wY, 1)
summary(glm(wY ~ t(wX) - 1, family="binomial"))
#ret = adagrad_sparse_sgd(X, y, rep(0, ncol(X)), rep(0, ncol(X)), 0.1)
