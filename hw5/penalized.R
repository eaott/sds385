########################################
# Penalized likelihood and soft thresholding
# (A)
########################################

########################################
# FIXME: run this for a few values.
########################################
par(mfrow=c(3,3))
yVals = c(8, 0, -1)
lambdaVals = c(0.5, 2, 10)
for (yVal in yVals) {
  for (lambdaVal in lambdaVals) {
    obj = list(y = yVal, lambda=lambdaVal)
    curve(0.5 * (obj$y - x)^2 + obj$lambda * abs(x), from=-15, to=15, xname="x",
          ylab=expression(paste("objective(", theta, ")")), xlab=expression(theta),
          main=bquote(y==.(obj$y) ~ ", " ~ lambda == .(obj$lambda)))
    S_lambda = sign(obj$y) * max(abs(obj$y) - obj$lambda, 0)
    abline(v=S_lambda)
    print(S_lambda)
    # hard threshold
    H_lambda = ifelse(obj$y >= obj$lambda, obj$y, 0)
    abline(v=H_lambda, col="red")
  }
}

########################################
# (B)
########################################
library(ggplot2)
library(reshape2)
n = 100
p = 0.1
thetas = 10 * (2 * rbinom(n, 1, 0.5) - 1) * rbinom(n, 1, p)
sigmas = rep(2, n)
zs = rnorm(n, thetas, sigmas)
lambdas = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 5, 10)
S_lambda = sapply(lambdas, function(lambda) {sign(zs) * pmax(abs(zs) - lambda * sigmas * sigmas, 0)})

data = data.frame(
  values = c(S_lambda),
  lmb = rep(lambdas, each=n),
  th = factor(rep(thetas, length(lambdas)))
)
par(mfrow=c(1,1))
ggplot(data, aes(x = th, y=values)) + geom_boxplot() +
  facet_wrap(~ lmb) +
  xlab(expression(theta)) + ylab(quote(widehat(theta)(y[i])))

mse = sapply(1:length(lambdas), function(i){a = S_lambda[, i] - thetas;
                              return(sum(a * a) / n)})
plot(mse ~ lambdas, type="l", log="xy")

