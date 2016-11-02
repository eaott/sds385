library(Matrix)
library(microbenchmark)
source("laplacian_methods.R")
source("graph_fused_lasso_methods.R")

fmri = read.csv("fmri_z.csv")
colnames(fmri) = NULL
fmri = as.matrix(fmri)
c.fmri = c(fmri)

D = makeD2_sparse(nrow(fmri), ncol(fmri))
L = t(D) %*% D

lambda = 2
C = Diagonal(128 ** 2) + lambda * L
n.iter = 200
# Direct solver for comparison
x = lu.solve(C, c.fmri)
admm.x = graphfusedlasso.solve(D, c.fmri, lambda = 2, n.iter = 100)

xMat = matrix(x, nrow = 128)
admmxMat = matrix(admm.x, nrow = 128)


par(mfrow = c(2, 2))
# cols = gray((0:50 / 50) ^ 0.67)
cols = colorspace::diverge_hcl(20)
image(fmri, col = cols, main = "raw")
image(xMat, col = cols, main = "lu")
image(admmxMat, col = cols, main = "admm")
sum(abs(admm.x) < 1e-15)

