library(Matrix)
library(ggplot2)

########################################
# Input
########################################
# Read in data and scale it.
x = read.csv("diabetesX.csv", header = TRUE)
y = read.csv("diabetesY.csv", header = FALSE)
colnames(x) = NULL
x = scale(x)
y = as.numeric(scale(y))

########################################
# Utility functions
########################################
# Numerical tolerance.
TOL = 1e-18

objective.fn = function(x, y, beta, lambda) {
  # Returns the lasso objective function value.
  #
  # Args:
  #   x: design matrix
  #   y: response vector
  #   beta: fitted parameters
  #   lambda: penalty multiplier
  return(sum((y - x %*% beta) ^ 2) + lambda * sum(abs(beta)))
}

is.converged = function(b1, b2, tol=1e-8) {
  # Simple convergence diagnostic for determining if b1 and b2 are
  # close enough to be converged.
  #
  # Args:
  #   b1: first vector
  #   b2: second vector
  #   tol: numerical tolerance
  #
  # Returns:
  #   Determines if the sum of squared differences between vectors is less
  #   then the tolerance.
  return(sum((b1 - b2) ^ 2) < tol)
}

########################################
# Algorithms
########################################
prox = function(x, y, n.iter, lambda.lasso, gamma) {
  # (Non-Accelerated) Proximal Gradient implementation.
  #
  # Args:
  #   x: design matrix
  #   y: response vector
  #   n.iter: maximum number of iterations
  #   lambda.lasso: penalty parameter for lasso
  #   gamma: tuning parameter; somewhat like a step multiplier
  #
  # Returns:
  #   Returns a list with elements:
  #   $n: number of iterations run
  #   $beta: final parameters of the fit
  #   $objective: value of the objective function for each iteration.

  # Scale up lambda for this algorithm.
  lambda = 2 * length(y) * lambda.lasso

  # Cache X^\top X and X^\top y
  xtx = crossprod(x)
  xy = t(x) %*% y

  # Variables for algorithm and return values.
  beta = rep(0, ncol(x))
  objective = rep(0, n.iter)

  # Use for convergence, determining how many data points to return.
  n.end = n.iter

  for (t in 1:n.iter) {
    # Store to check convergence at end of loop.
    oldbeta = beta

    # Main algorithm
    u = beta - 2 * gamma * (xtx %*% beta - xy)
    beta = sign(u) * pmax(abs(u) - lambda * gamma, 0)

    # Store objective value, check convergence.
    objective[t] = objective.fn(x, y, beta, lambda)
    if (is.converged(beta, oldbeta, TOL)) {
      n.end = t
      break
    }
  }
  return(list(n = n.end,
              beta = beta,
              objective = objective[1:n.end]))
}

accel.prox = function(x, y, n.iter, lambda.lasso, gamma) {
  # Accelerated Proximal Gradient implementation.
  #
  # Args:
  #   x: design matrix
  #   y: response vector
  #   n.iter: maximum number of iterations
  #   lambda.lasso: penalty parameter for lasso
  #   gamma: tuning parameter; somewhat like a step multiplier
  #
  # Returns:
  #   Returns a list with elements:
  #   $n: number of iterations run
  #   $beta: final parameters of the fit
  #   $objective: value of the objective function for each iteration.

  # Scale up lambda for this algorithm.
  lambda = 2 * length(y) * lambda.lasso

  # Cache X^\top X and X^\top y
  xtx = crossprod(x)
  xy = t(x) %*% y

  # Variables for algorithm and return values.
  beta.acc = rep(0, ncol(x))
  objective.acc = rep(0, n.iter)
  z.acc = beta.acc

  # The "momentum" value
  st = 1

  # Use for convergence, determining how many data points to return.
  n.end = n.iter

  for (t in 1:n.iter) {
    # Main algorithm
    oldbeta = beta.acc
    u.acc = z.acc - 2 * gamma * (xtx %*% z.acc - xy)
    beta.acc = sign(u.acc) * pmax(abs(u.acc) - lambda * gamma, 0)
    st1 = (1 + (1 + 4 * st ^ 2) ^ 0.5) / 2
    z.acc = beta.acc + (st - 1) / st1 * (beta.acc - oldbeta)
    st = st1

    # Store objective value, check convergence.
    objective.acc[t] = objective.fn(x, y, beta.acc, lambda)
    if (is.converged(beta.acc, oldbeta, TOL)) {
      n.end = t
      break
    }
  }
  return(list(n = n.end,
              beta = beta.acc,
              objective = objective.acc[1:n.end]))
}

admm = function(x, y, n.iter, lambda.lasso, rho) {
  # ADMM implementation
  #
  # Args:
  #   x: design matrix
  #   y: response vector
  #   n.iter: maximum number of iterations
  #   lambda.lasso: penalty parameter for lasso
  #   rho: tuning parameter; somewhat like a step multiplier
  #
  # Returns:
  #   Returns a list with elements:
  #   $n: number of iterations run
  #   $beta: final parameters of the fit
  #   $objective: value of the objective function for each iteration.

  # Scale up lambda for this algorithm.
  lambda = length(y) * lambda.lasso

  # Cache X^\top X and X^\top y
  xtx = crossprod(x)
  xy = t(x) %*% y

  # Variables for algorithm and return values.
  beta.admm = rep(0, ncol(x))
  objective.admm = rep(0, n.iter)
  z.admm = beta.admm
  u.admm = beta.admm
  xtx.cache = xtx + rho * diag(nrow(xtx))

  # Use for convergence, determining how many data points to return.
  n.end = n.iter

  for (t in 1:n.iter) {
    # Store to check convergence at end of loop.
    oldbeta = beta.admm

    # Main algorithm
    beta.admm = solve(xtx.cache, xy + rho * (z.admm - u.admm))
    temp = beta.admm + u.admm
    z.admm = sign(temp) * pmax(abs(temp) - lambda / rho, 0)
    u.admm = u.admm + beta.admm - z.admm

    ########################################
    # FIXME: should use the ADMM-specific convergence
    # criteria instead.
    ########################################

    # Store objective value, check convergence.
    objective.admm[t] = objective.fn(x, y, beta.admm, 2 * lambda)
    if (is.converged(beta.admm, oldbeta, TOL)) {
      n.end = t
      break
    }
  }
  return(list(n = n.end,
              beta = beta.admm,
              objective = objective.admm[1:n.end]))
}


########################################
# Plotting
########################################

# Using C_p from glmnet's default behavior, get lambda \approx 0.03
# lambda = 0.01   # admm .... prox/accel.prox
# lambda = 0.03   # admm ... prox .. accel.prox
# lambda = 0.07   # admm . accel.prox . prox
# lambda = 0.1    # prox .. accel.prox ... admm
lambda = 0.03

# Get result from glmnet to compare against because we can.
glmResult = glmnet::glmnet(x, y, alpha = 1, lambda = lambda)

# sidebar: test a decreasing sequence of lambdas instead for comparison
glm.n = 100
glm.scale = 20
glm.exp = 0.5
lambda.seq = seq((glm.scale * lambda) ^ glm.exp,
                 lambda ^ glm.exp,
                 length.out = glm.n) ^ (1/glm.exp)
glmResult2 = glmnet::glmnet(x, y, alpha = 1, lambda = lambda.seq)

# Now, actually get results from algorithms (lambda scaled inside functions).
n.iter = 1000
rho = 100
gamma = 0.0001
admm.result = admm(x, y, n.iter, lambda, rho)
prox.result = prox(x, y, n.iter,  lambda, gamma)
accel.prox.result = accel.prox(x, y, n.iter,  lambda, gamma)

# Format data to make ggplot easy to display the objective function
data = data.frame(objective = c(admm.result$objective,
                                prox.result$objective,
                                accel.prox.result$objective),
                  # plot against time
                  iter = c(1:admm.result$n,
                           1:prox.result$n,
                           1:accel.prox.result$n),
                  # for vertical bars
                  n = c(rep(admm.result$n, admm.result$n),
                        rep(prox.result$n, prox.result$n),
                        rep(accel.prox.result$n, accel.prox.result$n)),
                  # grouping
                  algorithm = c(rep("admm", admm.result$n),
                                rep("prox", prox.result$n),
                                rep("accel.prox", accel.prox.result$n)))
ggplot(data) + geom_line(aes(iter, objective, group = algorithm, col = algorithm)) +
  geom_vline(aes(xintercept = n, col = algorithm, group = algorithm), show.legend = FALSE) +
  scale_x_log10()

# Now, look at the coefficients instead of the objective.
res = round(data.frame(as.numeric(glmResult$beta),
                       as.numeric(glmResult2$beta[, glm.n]),
                       admm.result$beta,
                       prox.result$beta,
                       accel.prox.result$beta), 8)
names(res) = c("glm", "glm v.2", "admm", "prox", "accel.prox")

# Plot the coefficients in the same image
par(mfrow = c(1, 3))
plot(res$prox ~ res$glm,
     pch = 17, col = "blue", xlab = "GLM", ylab = "proximal gradient",
     main = "prox vs. GLM", asp = 1)
abline(a = 0, b = 1, h = 0, v = 0)
plot(res$accel.prox ~ res$glm,
     pch = 18, col = "red", xlab = "GLM", ylab = "Accelerated proximal gradient",
     main = "accel vs. GLM", asp = 1)
abline(a = 0, b = 1, h = 0, v = 0)
plot(res$admm ~ res$glm,
     pch = 19, col = "green", xlab = "GLM", ylab = "ADMM",
     main = "ADMM vs. GLM", asp = 1)
abline(a = 0, b = 1, h = 0, v = 0)

# also show coefficients textually
res
head(res, 10)

