library(Matrix)
source("makeD2_sparse.R")

########################################
# FIXME: dummy
########################################

graphfusedlasso.solve = function(A, b, n.iter = 200, tol = 1e-20) {
  if (!("dgCMatrix" %in% class(A))) {
    stop("Design matrix A should be a dgCMatrix")
  }

  x = rep(0, length(b))
  for (i in 1:n.iter) {
    x.old = x

    ########################################
    # FIXME: need the ADMM convergence conditions
    ########################################

    if (sum(abs(x.old - x) ^ 2) < tol) {
      #print(paste("gs:", i, "iterations"))
      break
    }
  }
  return(x)
}
