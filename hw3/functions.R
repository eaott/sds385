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
  w = 1 / (1 + exp(-as.matrix(x) %*% as.matrix(b)))
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

is.sufficient.decrease = function(X, Y, beta, m,
                                  oldLogL, c, tA, gL, pk) {
  newLogL = -logL(X, beta + tA * pk, Y, m)
  linearTerm = c * tA * sum(gL * pk)
  oldCond = -oldLogL + linearTerm
  #print(paste(newLogL, linearTerm, -oldLogL, oldCond))
  return(newLogL <= oldCond)
}

line_search = function(X, Y, m = 1, iter = 100,
                       c = 0.5, rho = 0.5, alpha = 0.1) {
  P = ncol(X)
  beta = rep(0, P)  # Could use random starting conditions as well.
  logLike = rep(0, iter)
  beta1 = rep(0, iter)
  for (i in 1:iter) {
      gL = gradL(X, beta, Y, m)
      pk = -gL
      oldLogL = logL(X, beta, Y, m)
      tempAlpha = alpha
      while (!is.sufficient.decrease(X, Y, beta, m,
                                     oldLogL, c, tempAlpha,
                                     gL, pk)) {
        tempAlpha = tempAlpha * rho
      }
      beta = beta + tempAlpha * pk
      logLike[i] = oldLogL
      beta1[i] = beta[1]
  }
  return(list(nll = -logLike,
              beta = beta,
              beta1 = beta1))
}


########################################
# Quasi-Newton
########################################

quasi_newton = function(X, Y, m = 1, iter = 100,
                        c = 0.5, rho = 0.5, alpha = 0.1) {
  P = ncol(X)

  beta = rep(0, P)  # Could use random starting conditions as well.
  oldBeta = beta
  hessian.inv = diag(P)

  logLike = rep(0, iter)
  beta1 = rep(0, iter)

  gL = rep(0, P)
  oldGl = rep(0, P)
counter = 0
  for (i in 1:iter) {
    # Update the old and new gradients
    oldGl = gL
    gL = gradL(X, beta, Y, m)

    # Use the previous hessian to compute the direction
    pk = -hessian.inv %*% gL

    # compute the log-likelihood
    oldLogL = logL(X, beta, Y, m)
    tempAlpha = alpha
    while (!is.sufficient.decrease(X, Y, beta, m,
                                   oldLogL, c, tempAlpha,
                                   gL, pk)) {
      tempAlpha = tempAlpha * rho
    }
    print(tempAlpha)
    # update the betas
    oldBeta = beta
    beta = beta + tempAlpha * pk

    # Convenience variables to follow the notation in the textbook
    sk = beta - oldBeta
    yk = gL - oldGl

    rhok = 1 / sum(sk * yk)

    ########################################
    # FIXME: I tried this, and it definitely
    # did not work as intended. It
    # basically blew everything up.
    ########################################
    # Update based on first update, see equation (6.20) in textbook
    #if (i == 1) {
    #  hessian.inv = sum(yk * yk) / sum(yk * sk) * hessian.inv
    #}


    ########################################
    # FIXME: rhok seems to be blowing up...
    # This ends up leading to NaN being
    # everywhere in the hessian and then
    # making life hard
    #
    # In the textbook, it talks about this on page 143 (bottom)
    # stating that this form of backtracking line search doesn't
    # work with this hessian model (works with gradient because
    # pk = -gradL).
    #
    # It says the update can be skipped, but it's not recommended
    # and that instead, a damped BFGS should be used (chapter 18).
    #
    # For now, I'm going to try the skipped method when 1/rhok is
    # negative or close to zero.
    #
    # In 2000 iterations, this executed 18 times... which means
    # it's not getting hardly any information about the
    # curvature. However, it does seem to converge
    # a bit faster than the steepest descent line search.
    ########################################
    if (sum(sk * yk) > 1e-20) {
      counter = counter + 1
      hessian.inv = (diag(P) - rhok * sk %*% t(yk)) %*%
        hessian.inv %*%
        (diag(P) - rhok * yk %*% t(sk)) + rhok * sk %*% t(sk)
    }

    logLike[i] = oldLogL
    beta1[i] = beta[1]
  }
  print(counter)
  return(list(nll = -logLike,
              beta = beta,
              beta1 = beta1))
}



