library(Rcpp)
library(RcppEigen)
Rcpp::sourceCpp("example.cpp")
wdbc = read.csv("wdbc.csv")

# Switch between using a scaled and non-scaled design matrix
X = scale(wdbc[ , 3:12])
#X = wdbc[ , 3:12]
# Add a column of 1s.
X = unname(as.matrix(cbind(1, X)))

Y = wdbc[ , 2]
# Y: replace M with 1 and B with 0.
Y = (Y == "M") * 1
m = rep(1, nrow(X))

nIter = 10000
eta = 1.0
ret = adagrad_gd(X, Y, m, nIter, eta)
beta = ret$beta
nll = ret$logL

gold.standard.betas = coef(glm(Y ~ X - 1, family="binomial"))
logL(X, gold.standard.betas, Y, m)
logL(X, beta, Y, m)


########################################
# SGD variant
########################################
ret = adagrad_sgd(X, Y, m, rep(0, ncol(X)), rep(0, ncol(X)), 0.1)
ll = ret$logL
for (i in 1:200) {
  ret = adagrad_sgd(X, Y, m, ret$beta,ret$adaDiag, 0.1)
  ll = c(ll, ret$logL)
}
ret$beta

fll = ll
fll[1] = nrow(X) * fll[1]
for (i in 2:length(ll)) {
  # fll[i] = 0.99 * fll[i - 1] + 0.01 * nrow(X) * fll[i]
  fll[i] = ((i - 1) * fll[i - 1] + nrow(X) * fll[i]) / i
}

plot(fll, type="l")
