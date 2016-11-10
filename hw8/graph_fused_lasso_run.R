library(colorspace)
library(grid)
library(gridExtra)
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
########################################
# FIXME: why is admm not giving exact zero values?
#
# Um, so James' code has a final line
# x[y == 0] = 0 to induce the sparsity...
########################################

admm = graphfusedlasso.solve(D, c.fmri, lambda = lambda,
                               n.iter = n.iter,
                               tol.abs = 1e-8,
                               tol.rel = 1e-5)

xMat = matrix(x, nrow = 128)
gsxMat = matrix(gs.x, nrow = 128)
admmxMat = matrix(admm$x, nrow = 128)

sparseLocations = fmri != 0

par(mfrow = c(2, 2))
layout(matrix(c(1,2,3,4), ncol = 2))
# cols = gray((0:50 / 50) ^ 0.67)
cols = diverge_hcl(40)

get_image = function(mat, name) {
  sub = ""
  xlab = ""
  ylab = ""
  scales = list(y = list(at = NULL), x = list(at = NULL))
  cuts = 15
  return(image(mat, main = name,
               sub = sub,
               xlab = xlab,
               ylab = ylab,
               scales = scales,
               cuts = cuts))
}

i1 = get_image(Matrix(fmri), "raw")
i2 = get_image(Matrix(xMat * sparseLocations), "lu")
i3 = get_image(Matrix(gsxMat * sparseLocations), "gauss")
i4 = get_image(Matrix(admmxMat * sparseLocations), "admm")
grid.arrange(i1, i2, i3, i4, ncol=2)



c(length(levels(factor(admm$r[ , 1]))), length(admm$r))
sum(abs(D %*% admm$x - admm$r) ^ 2)
admm$n
