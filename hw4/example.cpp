#include <RcppEigen.h>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::Matrix2d;
using Eigen::VectorXd;
using Eigen::MatrixXd;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
double TOL = 1e-6;


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

// [[Rcpp::export]]
Matrix2d example() {
  Matrix2d m;
  m(0, 0) = 3;
  m(1, 0) = 2.5;
  m(0, 1) = -1;
  m(1, 1) = m(1, 0) + m(0, 1);
  return m;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */

