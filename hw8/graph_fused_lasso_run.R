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
x = lu.solve(C, c(fmri))



