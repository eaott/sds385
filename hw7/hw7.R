library(Matrix)
library(ggplot2)

x = read.csv("diabetesX.csv", header = TRUE)
y = read.csv("diabetesY.csv", header = FALSE)
x = data.matrix(x, rownames.force = FALSE)
y = data.matrix(y, rownames.force = FALSE)

varnames = colnames(x)
colnames(x) = NULL



########################################
# ADMM (new)
########################################
admm = function(x, y, n.iter, lambda, rho) {
  xtx = crossprod(x)
  xy = t(x) %*% y

  beta.admm = rep(0, ncol(x))
  objective.admm = rep(0, n.iter)
  z.admm = beta.admm
  u.admm = beta.admm

  xtx.inv = solve(xtx + rho * diag(nrow = nrow(xtx)))

  for (t in 1:n.iter) {
    beta.admm = xtx.inv %*% (xy + rho * (z.admm - u.admm))
    temp = beta.admm + u.admm
    z.admm = sign(temp) * pmax(abs(temp) - lambda / rho, 0)
    u.admm = u.admm + beta.admm - z.admm
    objective.admm[t] = sum((y - x %*% beta.admm) ^ 2) + lambda * sum(abs(beta.admm))
  }
  return(list(beta=beta.admm,
              objective=objective.admm))
}


########################################
# Accelerated Proximal Gradient (last week)
########################################
accel.prox = function(x, y, n.iter, lambda, gamma) {
  xtx = crossprod(x)
  xy = t(x) %*% y

  beta.acc = rep(0, ncol(x))
  objective.acc = rep(0, n.iter)
  z.acc = beta.acc
  st = 1
  for (t in 1:n.iter) {
    u.acc = z.acc - 2 * gamma * (xtx %*% z.acc - xy)
    oldbeta = beta.acc
    beta.acc = sign(u.acc) * pmax(abs(u.acc) - lambda * gamma, 0)
    st1 = (1 + (1 + 4 * st ^ 2) ^ 0.5) / 2
    z.acc = beta.acc + (st - 1) / st1 * (beta.acc - oldbeta)
    st = st1
    objective.acc[t] = sum((y - x %*% beta.acc) ^ 2) + lambda * sum(abs(beta.acc))
  }
  return(list(beta=beta.acc,
              objective=objective.acc))
}

########################################
# Plotting
########################################
n.iter = 1000
rho = 10
lambda = 10
gamma = 0.001

admm.result = admm(x, y, n.iter, lambda, rho)
prox.result = accel.prox(x, y, n.iter, lambda, gamma)
data = data.frame(objective = c(admm.result$objective,
                                prox.result$objective),
                  iter = rep(1:n.iter, 2),
                  algorithm = rep(c("admm", "prox"), each = n.iter))

ggplot(data) + geom_line(aes(iter, objective, group = algorithm, col = algorithm))


