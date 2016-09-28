#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::Matrix2d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace Rcpp;

double TOL = 1e-10;
double EPSILON = 1e-8;

// [[Rcpp::export]]
VectorXd weightCpp(MatrixXd x, VectorXd beta) {
  // Create a column of ones for the quotient
  if (x.cols() != beta.rows()) {
    std::cout << "Wrong dimensions in weightCpp" << std::endl;
  }
  VectorXd expPart = (-1.0 * x * beta).array().exp();
  VectorXd denominator = 1.0 + expPart.array();
  VectorXd weights = denominator.cwiseInverse();
  return weights.cwiseMax(TOL).cwiseMin(1 - TOL);
}

// [[Rcpp::export]]
VectorXd weightVCpp(VectorXd x, VectorXd beta) {
  // Create a column of ones for the quotient
  if (x.rows() != beta.rows() || x.cols() != beta.cols()) {
    std::cout << "Wrong dimensions in weightVCpp" << std::endl;
  }
  VectorXd expPart = (-1.0 * x.transpose() * beta).array().exp();
  VectorXd denominator = 1.0 + expPart.array();
  VectorXd weights = denominator.cwiseInverse();
  return weights.cwiseMax(TOL).cwiseMin(1 - TOL);
}

// [[Rcpp::export]]
VectorXd gradLCpp(MatrixXd x, VectorXd b, VectorXd y, VectorXd m) {
  if (x.cols() != b.rows() || x.rows() != y.rows() || x.rows() != m.rows()) {
    std::cout << "Wrong dimensions somewhere" << std::endl;
  }
  VectorXd weightPart = y.array() - m.array() * weightCpp(x, b).array();

  // What I want is asdf[i] * x[i, ]
  for (int r = 0; r < x.rows(); r++) {
    x.row(r) *= weightPart[r];
  }
  return -x.colwise().sum();
}

/**
 * FIXME
 * Can add a double lambda parameter and add the L2 regularization
 */

// [[Rcpp::export]]
VectorXd gradWCpp(MatrixXd x, VectorXd b, VectorXd y, VectorXd m, VectorXd weights) {
  if (x.rows() != y.rows() || x.rows() != m.rows()
        || weights.rows() != x.rows()) {
    std::cout << "Wrong dimensions somewhere" << std::endl;
  }
  VectorXd weightPart = y.array() - m.array() * weights.array();

  // What I want is asdf[i] * x[i, ]
  for (int r = 0; r < x.rows(); r++) {
    x.row(r) *= weightPart[r];
  }
  VectorXd ret = -x.colwise().sum();
  //ret = ret + lambda * 2 * b;
  return ret;
}

// [[Rcpp::export]]
double logLikelihood(MatrixXd x, VectorXd b, VectorXd y, VectorXd m) {
  VectorXd p = weightCpp(x, b);
  double f = (y.array() * p.array().log() + (m - y).array() * (1.0 - p.array()).log()).sum();
  return f;
}

// [[Rcpp::export]]
List adagrad_gd(MatrixXd x, VectorXd y, VectorXd m, int nIter, double eta) {
  /* Proof of concept for the AdaGrad solution.
   *
   * In this instance, this does the full thing for the full design
   * matrix, rather than something stochastic.
   */

  const int N = x.rows();
  const int P = x.cols();
  if (y.rows() != N || m.rows() != N) {
    std::cout << "Wrong dimensions" << std::endl;
  }
  // Create matrices here.
  VectorXd beta = VectorXd::Zero(P, 1);
  VectorXd gradient = VectorXd::Zero(P, 1);
  VectorXd adaDiag = VectorXd::Zero(P, 1);
  VectorXd denominator = VectorXd::Zero(P, 1);
  VectorXd weights = VectorXd::Zero(P, 1);

  VectorXd logL = VectorXd::Zero(nIter, 1);
  for (int i = 0; i < nIter; i++) {
    weights = weightCpp(x, beta);
    // Compute gradient on full data set
    gradient = gradLCpp(x, beta, y, m);
    // Compute the squares of the elements of the gradient and use that for the hessian-ish-thing
    adaDiag = adaDiag + (gradient.array() * gradient.array()).matrix();
    // Add numerical tolerance, take 1/sqrt(A) to multiply
    denominator = (adaDiag.array() + EPSILON).matrix().cwiseSqrt().cwiseInverse();
    // Do the update
    beta = beta - eta * denominator.cwiseProduct(gradient);
    // Store the log-likelihood
    logL[i] = logLikelihood(x, beta, y, m);
  }


  return List::create(
    Named("beta") = beta,
    Named("logL") = logL
  );
}

/**
 * FIXME
 * Here's the deal with the following function. It doesn't maintain the
 * old value of the log-likelihood. Could just return the new value
 * every time and post-process, or could provide a default value
 * of 0 and re-load it.
 */

// [[Rcpp::export]]
List adagrad_sgd(MatrixXd x, VectorXd y, VectorXd m,
                 VectorXd beta, VectorXd adaDiag,
                 double eta) {
  /*
   * In this case, we need additional parameters to keep the running
   * values of beta and adaDiag going.
   *
   * Just go through the available values of x?
   */

  const int N = x.rows();
  const int P = x.cols();
  if (y.rows() != N || m.rows() != N || beta.rows() != P || adaDiag.rows() != P) {
    std::cout << "Wrong dimensions" << std::endl;
  }
  // Create matrices here.
  VectorXd gradient = VectorXd::Zero(P, 1);
  VectorXd denominator = VectorXd::Zero(P, 1);
  VectorXd weights = VectorXd::Zero(P, 1);

  VectorXd logL = VectorXd::Zero(N, 1);

  for (int i = 0; i < N; i++) {
    // Compute weights only once.
    weights = weightVCpp(x.row(i), beta);
    std::cout << weights << std::endl;
    // Compute gradient on full data set
    gradient = gradWCpp(x.row(i), beta, y.row(i), m.row(i), weights);
    // Compute the squares of the elements of the gradient and use that for the hessian-ish-thing
    adaDiag = adaDiag + (gradient.array() * gradient.array()).matrix();
    // Add numerical tolerance, take 1/sqrt(A) to multiply
    denominator = (adaDiag.array() + EPSILON).matrix().cwiseSqrt().cwiseInverse();
    // Do the update
    beta = beta - eta * denominator.cwiseProduct(gradient);

    // Store the log-likelihood
    logL[i] = logLikelihood(x.row(i), beta, y.row(i), m.row(i));
  }


  return List::create(
    Named("adaDiag") = adaDiag,
    Named("beta") = beta,
    Named("logL") = logL
  );
}
