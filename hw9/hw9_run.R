library(Matrix)
u = matrix(c(1, 1,
             0, 1,
             1, 0), ncol = 2, byrow = TRUE)

v = matrix(c(0, 1,
             0, 1,
             1, 0,
             1, 1,
             1, 0), ncol = 2, byrow = TRUE)

d = c(0.2, 5, 1)
X = d * u %*% t(v)
rankMatrix(X)
# X = matrix(c(1e-5, 1, 0,
#          0, 1, 0,
#          0, 1, 0), nrow = 3)

lambda = 1
n.iter = 2000
TOL = 1e-30
res = pmd(X, lambda, lambda, 2, n.iter = n.iter, TOL = TOL)
print(res)
image(rbind(X, X - res$resid), col = rainbow(20))
#print(cbind(X, res$u %*% t(res$v) * res$d))
resid = X - res$u %*% t(res$v) * res$d
a = pmd(resid, lambda, n.iter = n.iter)
print(resid - a$u %*% t(a$v) * a$d)


data = read.csv("social_marketing.csv", header = TRUE)
data = data[ , 2:ncol(data)]
for (i in 1:ncol(data)) {
  data[ , i] = data[ , i] / max(data[ , i])
}
kinds = names(data)
names(data) = NULL
data = Matrix::Matrix(data = as.matrix(data))

# For checking determinism up to random seed.
set.seed(8675309)

lambda1 = 1e+1
lambda2 = 1.3
rank = 30
n.iter = 20
TOL = 1e-50
res = pmd(data, lambda1, lambda2, rank, n.iter = n.iter, TOL = TOL)
sum(res$u[,1] == 0)

for (i in 1:rank) {
  print(paste(kinds[which(res$v[ , order(-res$d)][ , i] != 0)], collapse = " "))
}
image(abs(res$v) > 0)
