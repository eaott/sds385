library(Matrix)
library(ggplot2)

x = read.csv("diabetesX.csv", header = TRUE)
y = read.csv("diabetesY.csv", header = FALSE)
x = data.matrix(x, rownames.force = FALSE)
y = data.matrix(y, rownames.force = FALSE)

varnames = colnames(x)
colnames(x) = NULL
xtx = crossprod(x)
xy = t(x) %*% y

n.iter = 1000
gamma = 0.001
lambda = 10
beta = rep(0, ncol(x))
objective = rep(0, n.iter)
for (t in 1:n.iter) {
  u = beta - 2 * gamma * (xtx %*% beta - xy)
  beta = sign(u) * pmax(abs(u) - lambda * gamma, 0)
  objective[t] = sum((y - x %*% beta) ^ 2) + lambda * sum(abs(beta))
}



########################################
# Accelerated proximal gradient
########################################

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


data = data.frame(old = objective,
                  new = objective.acc,
                  x = 1:n.iter)
ggplot(data) + geom_line(aes(x, old)) +
  geom_line(aes(x, new), col = "red")


