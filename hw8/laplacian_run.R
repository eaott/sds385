library(Matrix)
library(microbenchmark)
library(RColorBrewer)
source("laplacian_methods.R")
fmri = read.csv("fmri_z.csv")
colnames(fmri) = NULL
fmri = as.matrix(fmri)
c.fmri = c(fmri)

D = makeD2_sparse(nrow(fmri), ncol(fmri))
L = t(D) %*% D

lambda = 2
C = Diagonal(128 ** 2) + lambda * L
n.iter = 200
x = lu.solve(C, c.fmri)
gs.x = gaussseidel.solve(C, c.fmri, n.iter = n.iter)
jac.x = jacobi.solve(C, c.fmri, n.iter = n.iter)
# Plot
xMat = matrix(x, nrow = 128)
gsxMat = matrix(gs.x, nrow=128)
jacxMat = matrix(jac.x, nrow=128)

par(mfrow = c(2, 2))
# cols = gray((0:50 / 50) ^ 0.67)
cols = colorspace::diverge_hcl(20)
image(fmri, col = cols, main = "raw")
image(xMat, col = cols, main = "lu")
image(gsxMat, col = cols, main = "gs")
image(jacxMat, col = cols, main = "jac")

########################################
# FIXME: Use the sparse matrix image function instead.
########################################


# HOBBITSES ARE TRICKSEY
# This is a completely unfair comparison as-is.
# the LU decomposition actually gets stored back into C (a property
# of the Matrix package). After lu(C) has been called,
# C@factors$LU caches the decomposition.
# I found this out after lu.solve was finishing in 8 ms and getting the right
# answer.
microbenchmark::microbenchmark(lu.solve(C, c.fmri),
                               gaussseidel.solve(C, c.fmri, n.iter = n.iter),
                               jacobi.solve(C, c.fmri, n.iter = n.iter),
                               times = 10L)

# Added a flag to lu.solve to handle this with correct scoping
microbenchmark::microbenchmark(lu.solve(C, c.fmri, useCache = FALSE),
                               gaussseidel.solve(C, c.fmri, n.iter = n.iter),
                               jacobi.solve(C, c.fmri, n.iter = n.iter),
                               times = 10L)


