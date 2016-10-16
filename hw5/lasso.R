library(glmnet)
library(ggplot2)
library(grid)
library(gridExtra)
library("scales")
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

x = read.csv("diabetesX.csv", header = TRUE)
y = read.csv("diabetesY.csv", header = FALSE)
x = data.matrix(x, rownames.force = FALSE)
y = data.matrix(y, rownames.force = FALSE)
fit = glmnet(x, y)

dim(fit$beta)



# 64 x 100
plotdata = t(fit$beta)
plotdata = as.matrix(plotdata)
rownames(plotdata) = NULL
repLambdas = rep(fit$lambda, nrow(fit$beta))
plotdata = cbind(melt(plotdata), repLambdas)

b = y - x %*% fit$beta
base.mse = colMeans(b * b)

n = nrow(x)
ind = sample(1:n, n)
folds = 26 # n / 26 is an integer
mses = c()
########################################
# FIXME: use the `cut` function instead for indices
########################################

for (f in 1:folds) {
  starti = n / folds * (f - 1) + 1
  endi = n / folds * f
  test = ind[(starti:endi)]
  train = ind[-(starti:endi)]
  xtrain = x[train, ]
  xtest = x[test, ]
  ytrain = y[train, ]
  ytest = y[test, ]
  fit2 = glmnet(xtrain, ytrain, lambda=fit$lambda)

  b2 = ytest - xtest %*% fit2$beta
  mse = colMeans(b2 * b2)
  mses = rbind(mses, mse)
}
mses.mean = colMeans(mses)
mses.se = sapply(1:100, function(i){sd(mses[, i])})
# Use the global lambdas.
Cp = base.mse + 2 * fit$df / n * sapply(1:100, function(i){var(b[,i])})
myplot = data.frame(mse = mses.mean,
                    err = mses.se,
                    lmb = fit$lambda,
                    base = base.mse,
                    cp = Cp)


########################################
# FIXME: to get geom_vline to work nicely, need to add a column
# to myplot2 that basically transposes the type variable to the
# log(lambda[which.min(type)]) sort of thing.
########################################

mse_lambda = myplot$lmb[which.min(base.mse)]
moose_lambda = myplot$lmb[which.min(mses.mean)]
cp_lambda = myplot$lmb[which.min(Cp)]

myplot2 = data.frame(val = c(base.mse, mses.mean, Cp),
                     lmb = rep(fit$lambda, 3),
                     type = rep(c("mse", "moose", "cp"), each=length(fit$lambda)),
                     xint = rep(c(mse_lambda, moose_lambda, cp_lambda), each=length(fit$lambda)))
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
cols = ggplotColours(n = 3)
g2 = ggplot(myplot2, aes(x=log(lmb), y=val, col=type)) +
  geom_line() +
  xlim(c(max(log(fit$lambda)), min(log(fit$lambda)))) +
  geom_vline(xintercept=log(mse_lambda), col=cols[3]) +
  geom_vline(xintercept=log(moose_lambda), col=cols[2]) +
  geom_vline(xintercept=log(cp_lambda), col=cols[1]) +
  # geom_vline(aes(xintercept=log(xint), col=cols[as.factor(type)]), show.legend = FALSE) +
  xlab(expression(log(lambda))) +
  ylab("error") + theme(legend.position = c(.85, .65), legend.title=element_blank())

g1 = ggplot(plotdata, aes(x=log(repLambdas), y=value, group=Var2, col=Var2)) +
  geom_line() +
  geom_vline(xintercept=log(mse_lambda), col=cols[3]) +
  geom_vline(xintercept=log(moose_lambda), col=cols[2]) +
  geom_vline(xintercept=log(cp_lambda), col=cols[1]) +
  guides(group=FALSE, colour=FALSE) +
  xlab(expression(log(lambda))) + ylab(quote(widehat(beta)[lambda])) +
  xlim(c(max(log(fit$lambda)), min(log(fit$lambda))))


grid.arrange(g2, g1, ncol=1, heights=c(0.5, 1.5))
