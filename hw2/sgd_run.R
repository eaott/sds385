source("sgd_functions.R")
wdbc = read.csv("wdbc.csv")

# Switch between using a scaled and non-scaled design matrix
X = scale(wdbc[ , 3:12])
#X = wdbc[ , 3:12]
# Add a column of 1s.
X = unname(as.matrix(cbind(1, X)))

Y = wdbc[ , 2]
# Y: replace M with 1 and B with 0.
Y = as.matrix((Y == "M") * 1)

sgd(X, Y, 1, iter = 40000, step = 0.01, alpha = 0.5 / nrow(X))

########################################
# Part D still in progress
########################################
sgd_robbinsMonro(X, Y, 1, iter = 10000,
                 C = 1, t0 = 1, a = 0.55,
                 alpha = 0.5/nrow(X))
