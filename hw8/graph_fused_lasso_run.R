library(colorspace)
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
C = Diagonal(nrow(fmri) * ncol(fmri)) + lambda * L
n.iter = 600
# Direct solver for comparison
x = lu.solve(C, c.fmri)
gs.x = gaussseidel.solve(C, c.fmri, n.iter = n.iter)
admm = graphfusedlasso.solve(D, c.fmri, lambda = lambda,
                               n.iter = n.iter,
                               tol.abs = 1e-8,
                               tol.rel = 1e-5)

xMat = matrix(x, nrow = 128)
gsxMat = matrix(gs.x, nrow = 128)
admmxMat = matrix(admm$x, nrow = 128)


par(mfrow = c(2, 2))
# cols = gray((0:50 / 50) ^ 0.67)
cols = diverge_hcl(40)
image(fmri, col = cols, main = "raw")
image(xMat, col = cols, main = "lu")
image(gsxMat, col = cols, main = "gauss")
image(admmxMat, col = cols, main = "admm")

c(length(levels(factor(admm$r[ , 1]))), length(admm$r))
sum(abs(D %*% admm$x - admm$r) ^ 2)
admm$n
