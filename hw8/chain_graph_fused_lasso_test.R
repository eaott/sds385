library(ggplot2)
library(Matrix)
library(microbenchmark)
library(reshape2)
source("laplacian_methods.R")
source("graph_fused_lasso_methods.R")

N = 200
noise = rnorm(N, sd = 0.25)
orig = rep(0, N)
orig[50:80] = orig[50:80] + 5
orig[130:150] = orig[130:150] + 2
data = as.matrix(noise + orig)

D = makeD2_sparse(1, N)

L = t(D) %*% D

lambda = 10
C = Diagonal(nrow(data) * ncol(data)) + lambda * L
n.iter = 200


# Direct solver for comparison
lu.x = lu.solve(C, c(data))
admm.x = graphfusedlasso.solve(D, c(data), lambda = lambda, n.iter = 1000,
                               tol.abs = 1e-3, tol.rel = 1e-8)

plot_data = data.frame(
  iter = 1:N,
  orig = orig,
  data = data,
  lu = lu.x,
  admm = admm.x$z
)

melted_data = melt(plot_data, id.vars = c("iter", "data"))

ggplot(melted_data) + geom_point(aes(iter, data)) +
  geom_line(aes(iter, value, group = variable, color = variable))

length(unique(factor(admm.x$r[ , 1])))
unique(factor(admm.x$r[ , 1]))
