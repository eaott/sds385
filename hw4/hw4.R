library(ggplot2)
library(Rcpp)
library(RcppEigen)
Rcpp::sourceCpp("example.cpp")

exX = matrix(c(1,2,3,4,5,6,7,8,9), nrow=3)
exB = c(0.1, 0, 0.2)
exY = c(1,2,3)
exM = c(1,1,1)
# adagrad(exX, exY, exM, 10, 1.0)


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
                                  oldLogL, cc, tA, gL, pk) {
  newLogL = -logL(X, beta + tA * pk, Y, m)
  linearTerm = cc * tA * sum(gL * pk)
  oldCond = -oldLogL + linearTerm
  #print(paste(newLogL, linearTerm, -oldLogL, oldCond))
  return(newLogL <= oldCond)
}


########################################
# New functions (other from hw 1)
########################################
sgd = function(X, Y, m = 1, iter = 1000, step = 0.01, alpha=0.9, detail=F) {
  # Function that does the entire stochastic gradient descent
  # with a fixed step size and number of iterations. Trains on
  # only one data point at a time.
  #
  # Args:
  #   X: design matrix, N x P
  #   Y: response vector as a matrix, N x 1
  #   m: maximum value for i-th component of Y (can be scalar as a fixed value
  #      for all components of Y)
  #   iter: Number of iterations to run
  #   step: the fixed step size
  #   alpha: constant in the exponentially-weighted moving average
  #   detail: if TRUE, include the actual contributions for each data point
  #           in the resultant plot
  #
  # Returns:
  #   Plot of the negative log-likelihood as a function of iteration number,
  #   on a log-log scale. In black, the full negative log-likelihood. In green,
  #   the running average computed on each data point (on full data set scale).
  #   In red, the exponentially-weighted moving average.
  N = nrow(X)
  P = ncol(X)
  m = rep(m, N)
  beta = rep(0, P)  # Could use random starting conditions as well.
  logLike = rep(0, iter)
  windowLogLike = rep(0, iter)
  weightedLogLike = rep(-300, iter)
  contributions = rep(0, iter)

  getNewStep = TRUE
  currentStep = step

  for (i in 1:iter) {
    if (getNewStep) {
      ########################################
      # FIXME: something like this
      ########################################
      currentStep = step
      # TODO: Should sample only a mini-batch
      xCalib = X
      yCalib = Y
      mCalib = m
      oldLogL = logL(xCalib, beta, yCalib, mCalib)

      avgGrad = gradL(xCalib, beta, yCalib, mCalib) / nrow(xCalib)

      while (!is.sufficient.decrease(xCalib, yCalib, beta, mCalib,
                                   oldLogL, cc, currentStep,
                                   avgGrad, -avgGrad)) {
        currentStep = currentStep * rho
      }

    }



    pt = sample(1:N, 1)
    xSample = matrix(X[pt,], nrow = 1)
    ySample = matrix(Y[pt,], nrow = 1)
    mSample = m[pt]
    gL = gradL(xSample, beta, ySample, mSample)
    beta = beta - step * gL







    # Plot the actual log-likelihood, not just the sample.
    logLike[i] = logL(X, beta, Y, 1)

    # Log-likelihood samples
    contribution = N * logL(xSample, beta, ySample, 1)
    contributions[i] = contribution
    windowLogLike[i] = contribution
    # Simple moving average (scale previous entry to the sum of all previous
    # entries, add new entry, divide by new size)
    if (i > 1) {
      windowLogLike[i] = ((i - 1) * windowLogLike[i - 1] + windowLogLike[i]) / i
    }
    # For the exponentially-weighted result, start after 10% of iterations.
    # At this point, back-fill all the results to the current contribution
    # (a terrible guess, but if we ignore the burn-in, it's no different than
    # starting at the first value)
    if (i == round(iter / 10)) {
      weightedLogLike[1:i] = contribution
    }
    # After the burn-in, take a weighted average between the new contribution
    # and old one.
    if (i > round(iter / 10)) {
      weightedLogLike[i] = alpha * contribution +
        (1 - alpha) * weightedLogLike[i - 1]
    }
  }
  print(paste("Final value of full log-likelihood:", -logLike[iter]))
}
