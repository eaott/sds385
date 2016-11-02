library(Matrix)
library(microbenchmark)
source("laplacian_methods.R")
source("graph_fused_lasso_methods.R")


orig = c(1, 1, 1, 1, 1,
        2, 2, 1, 1, 3,
        2, 2, 1, 3, 3,
        5, 4, 4, 3, 4,
        5, 5, 4, 4 ,4)
data = orig + rnorm(length(orig), sd = 0.1)
data = matrix(data, nrow = 5)
D = makeD2_sparse(nrow(data), ncol(data))

L = t(D) %*% D

lambda = 2
C = Diagonal(nrow(data) * ncol(data)) + lambda * L
n.iter = 200
# Direct solver for comparison
x = lu.solve(C, c(data))
admm.x = graphfusedlasso.solve(D, c(data), lambda = lambda / length(orig), n.iter = 1000,
                               tol.abs = 1e-3, tol.rel = 1e-8)

par(mfrow = c(2, 2))
# cols = gray((0:50 / 50) ^ 0.67)
cols = colorspace::diverge_hcl(25)
image(matrix(orig, nrow = 5), col = cols, main = "orig")
image(matrix(data, nrow = 5), col = cols, main = "data")
image(matrix(x, nrow = 5), col = cols, main = "lu")
image(matrix(admm.x$x, nrow = 5), col = cols, main = "admm")
length(unique(factor(admm.x$r[ , 1])))
sum(abs(D %*% admm.x$x - admm.x$r) ^ 2)
admm.x$n
