source("functions.R")
wdbc = read.csv("wdbc.csv")

# Switch between using a scaled and non-scaled design matrix
X = scale(wdbc[ , 3:12])
#X = wdbc[ , 3:12]
# Add a column of 1s.
X = unname(as.matrix(cbind(1, X)))

Y = wdbc[ , 2]
# Y: replace M with 1 and B with 0.
Y = as.matrix((Y == "M") * 1)

P = ncol(X)
gold.standard.betas = as.numeric(coef(glm(Y ~ X - 1, family = "binomial")))
gold.standard = -logL(X, gold.standard.betas, Y, 1)

iter = 20000
ls = line_search(X, Y, m = 1, iter = iter, cc = 0.9, rho = 0.5, alpha = 0.1)
plot(ls$nll)
print(paste("ideal:", gold.standard, "mine:", tail(ls$nll, 1)))
print(gold.standard.betas)
print(ls$beta)
plot(ls$beta ~ gold.standard.betas, asp = 1)
plot(ls$beta1, type = "l")

########################################
# Quasi-Newton
########################################

qn = quasi_newton(X, Y, m = 1, iter = 200, cc = 0.1, rho = 0.1, alpha = 1)
plot(qn$nll)

