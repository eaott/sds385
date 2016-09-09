weight = function(x, b) {
  # Function to compute the weight for logistic regression.
  #
  # Args:
  #   x: Design matrix, N x P
  #   b: Regression parameters, P x 1
  #
  # Returns:
  #   The weight. If x is a matrix, will have same number of rows. Resultant
  #   weights will be in the range [1e-6, 1-1e-6] for stability
  w = (1/(1+exp(-as.matrix(x) %*% as.matrix(b))))
  w = pmax(w, 1e-6)
  w = pmin(w, 1 - 1e-6)
  return(w)
}


gradL = function(x, b, y, m) {
  # Computes the gradient of the negative log-likelihood with respect to the
  # regression parameters.
  #
  # Args:
  #   x: Design matrix, N x P
  #   b: Regression parameters, P x 1
  #   y: Response vector, N x 1. Each value y_i must be in {0, ..., m_i}
  #   m: max value vector, N x 1 (or a scalar) for use above
  #
  # Returns:
  #   The gradient as a P x 1 vector.

  # This will multiply each element into the corresponding row of x
  terms = as.numeric(y - m * weight(x, b)) * x
  # now, take the multiplied rows and sum them up creating
  # a 1xP vector that is just a c(...) so Px1
  gL = colSums(terms)
  return(-gL)
}

logL = function(x, b, y, m) {
  # Computes the log-likelihood.
  #
  # Args:
  #   x: Design matrix, N x P
  #   b: Regression parameters, P x 1
  #   y: Response vector, N x 1. Each value y_i must be in {0, ..., m_i}
  #   m: max value vector, N x 1 (or a scalar) for use above
  #
  # Returns:
  #   The log-likelihood.

  # Leverages the dbinom function which can take in vectors for each parameter
  # and do the computations accordingly, and cast the results into log form.
  return(sum(dbinom(y, m, weight(x, b), log = TRUE)))
}

########################################
# Part D
########################################

gradLNewton = function(x, b, y, m) {
  # Computes the gradient of the negative log-likelihood using Newton's method
  #
  # Args:
  #   x: Design matrix, N x P
  #   b: Regression parameters, P x 1
  #   y: Response vector, N x 1. Each value y_i must be in {0, ..., m_i}
  #   m: max value vector, N x 1 (or a scalar) for use above
  #
  # Returns:
  #   The gradient as a P x 1 vector.
  w = c(weight(x, b))  # Convert to vector instead of matrix
  hessian = t(x) %*% (m * w * (1-w) * x)
  gL = as.matrix(gradL(x, b, y, m))
  return(-solve(hessian, gL))  # Shortcut that uses a decomposition
}



