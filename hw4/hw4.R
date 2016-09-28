library(ggplot2)
library(Rcpp)
library(RcppEigen)
Rcpp::sourceCpp("example.cpp")

exX = matrix(c(1,2,3,4,5,6,7,8,9), nrow=3)
exB = c(0.1, 0, 0.2)
exY = c(1,2,3)
exM = c(1,1,1)
# adagrad(exX, exY, exM, 10, 1.0)





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

      avgGrad = gradLCpp(xCalib, beta, yCalib, mCalib) / nrow(xCalib)

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
    gL = gradLCpp(xSample, beta, ySample, mSample)
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
