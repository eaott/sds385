source("glm.R")
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
# newton = TRUE <==> iter ~ 100 is completely sufficient (even iter ~ 10)
# newton = FALSE <==> iter ~ 20000+
newton = TRUE
iter = 10

beta = rep(0, P)  # Could use random starting conditions as well.
logLike = rep(0, iter)
for (i in 1:iter) {
  if (newton) {
    ########################################
    # Part D
    ########################################
    beta = beta + gradLNewton(X, beta, Y, 1)
  } else {
    ########################################
    # Part B
    ########################################
    gL = gradL(X, beta, Y, 1)
    # For an X matrix that has been scaled
    beta = beta - 0.001 * gL
    # For an X matrix in the original scale
    # beta = beta - 0.001 * gL / sqrt(sum(gL ^ 2))
  }
  logLike[i] = logL(X, beta, Y, 1)
}
plot(logLike, type = "l")

# Particular results
# Scaled X, gradient
# 20000 iters
# LogL = -74.10653
#  [1] -0.2504395 -0.5149440  1.6259505 -0.9047349  5.4986217
# [6]  1.0729947 -0.5458339  0.9790104  2.3175757  0.4692481
# [11] -0.2735361

# Scaled X, newton
# 200 iters
# LogL = -73.06523
# [1]  0.46964720 -7.22158870  1.64995454 -1.73632394 14.00606994
# [6]  1.07358670 -0.07657184  0.67150243  2.58049248  0.44471461
# [11] -0.48075207
