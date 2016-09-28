#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::Matrix2d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::RowMajor;
using namespace Rcpp;
typedef SparseVector<double> SVd;
typedef SparseMatrix<double, RowMajor> SMDR;
typedef SparseMatrix<double> MyMatrix;

/*
 * This file is almost identical to the other one.
 * Differences include setting m=1 everywhere for the LR
 * in order to not use an extra array. It also handles the sparsity.
 */

double TOL = 1e-6;
//double EPSILON = 1e-8;

// [[Rcpp::export]]
double weightVCpp(VectorXd x, VectorXd beta) {
  // Create a column of ones for the quotient
  // std::cout << x.rows() << " " << x.cols() << " " << beta.rows() << " " << beta.cols() << std::endl;
  // if (x.rows() != beta.rows() || x.cols() != beta.cols()) {
  //   std::cout << "Wrong dimensions in weightVCpp" << std::endl;
  // }
  //VectorXd expPart = (-1.0 * x.transpose() * beta).array().exp();
  double expPart = exp(-1.0 * x.cwiseProduct(beta).sum());
  double denominator = 1.0 + expPart;
  double weight = 1.0 / denominator;
  return weight < TOL ? TOL : weight > 1 - TOL ? 1 - TOL : weight;
}


/**
* FIXME
* Can add a double lambda parameter and add the L2 regularization
* would need to convert to non-sparse because of beta.
*/

// [[Rcpp::export]]
VectorXd gradCpp(VectorXd x, VectorXd b, double y, double weight) {
  // The weight allows it to maybe be computed only once, rather than twice.
  // if (x.cols() != b.cols() || x.rows() != b.rows()) {
  //   std::cout << "Wrong dimensions somewhere" << std::endl;
  // }
  double weightPart = y - weight;
  VectorXd newX = x * weightPart;
  VectorXd ret = -x;
  //ret = ret + lambda * 2 * b;
  return ret;
}

// [[Rcpp::export]]
VectorXd gradSMCpp(VectorXd grad, SMDR x, int i, VectorXd b, double y, double weight) {
  // The weight allows it to maybe be computed only once, rather than twice.
  // if (x.cols() != b.cols() || x.rows() != b.rows()) {
  //   std::cout << "Wrong dimensions somewhere" << std::endl;
  // }
  double weightPart = y - weight;

  for (SMDR::InnerIterator it(x,i); it; ++it)
  {
    grad(it.col()) = -it.value() * weightPart;
  }
  //ret = ret + lambda * 2 * b;
  return grad;
}

// [[Rcpp::export]]
VectorXd reset(VectorXd z, SMDR x, int i) {
  // The weight allows it to maybe be computed only once, rather than twice.
  // if (x.cols() != b.cols() || x.rows() != b.rows()) {
  //   std::cout << "Wrong dimensions somewhere" << std::endl;
  // }

  for (SMDR::InnerIterator it(x,i); it; ++it)
  {
    z(it.col()) = 0;
  }
  //ret = ret + lambda * 2 * b;
  return z;
}


// [[Rcpp::export]]
double logLikelihood(VectorXd b, double y, double p) {
  double f = y * log(p) + (1 - y) * log(1 - p);
  return f;
}


/**
* FIXME
* Here's the deal with the following function. It doesn't maintain the
* old value of the log-likelihood. Could just return the new value
* every time and post-process, or could provide a default value
* of 0 and re-load it.
*/

// [[Rcpp::export]]
double test(MyMatrix xCWISE) {
  /*
  * In this case, we need additional parameters to keep the running
  * values of beta and adaDiag going.
  *
  * Just go through the available values of x?
  */
  // SMDR x(xCWISE);
  // std::cout << "finished restructuring" << std::endl;
  std::cout << "got to c++" << std::endl;

  const int N = xCWISE.cols(); // xCWISE is a P x N matrix!!!
  const int P = xCWISE.rows(); // xCWISE is a P x N matrix!!!
  VectorXd beta = VectorXd::Zero(P, 1);


  for (int k=0; k< N; ++k) {
    ///////////////////////////////////////////
    // this section computes the weight for a single sample.
    // it's set up this way because calling a function
    // was tested (including removing the Rcpp::export
    // and including adding the inline to the function definition)
    // and doing it all properly inline was at least 10x faster.
    // I honestly didn't wait long enough to see if it was faster
    // than that, and instead just killed it. Literally
    // running the same code in a function slowed it down like crazy.
    ///////////////////////////////////////////
    double sum = xCWISE.col(k).dot(beta);
    double expPart = exp(-1.0 * sum);
    double denominator = 1.0 + expPart;
    double weight = 1.0 / denominator;
    // update weight to be in [TOL, 1 - TOL].
    weight = (weight < TOL ? TOL : weight > 1 - TOL ? 1 - TOL : weight);
  }
  return 0.0;
}









// [[Rcpp::export]]
List adagrad_sparse_sgd(SparseMatrix<double> xCWISE,
                        VectorXd y,
                        VectorXd beta,
                        VectorXd adaDiag,
                        double eta) {
  /*
  * In this case, we need additional parameters to keep the running
  * values of beta and adaDiag going.
  *
  * Just go through the available values of x?
  */
  SMDR x(xCWISE);
  std::cout << "finished restructuring" << std::endl;
  const int N = x.rows();
  const int P = x.cols();
  if (y.rows() != N || beta.rows() != P || adaDiag.rows() != P) {
    std::cout << "Wrong dimensions" << std::endl;
  }
  // Create matrices here.
  VectorXd zeros = VectorXd::Zero(P, 1);
  VectorXd gradient = VectorXd::Zero(P, 1);
  VectorXd denominator = VectorXd::Zero(P, 1);
  double weight = 0;

  VectorXd logL = VectorXd::Zero(N, 1);
  //
  //   for (int k=0; k<x.outerSize(); ++k) {
  //     for (SparseMatrix<double>::InnerIterator it(x,k); it; ++it)
  //     {
  //       it.value();
  //       it.row();   // row index
  //       it.col();   // col index (here it is equal to k)
  //       it.index(); // inner index, here it is equal to it.row()
  //     }
  //   }


  for (int i = 0; i < 1; i++) {


    // Old way that kinda works if I use MatrixXd x.
    // VectorXd row = x.row(i);
    // weight = weightVCpp(row, beta);
    // gradient = gradCpp(row, beta, y(i), weight);
    // adaDiag = adaDiag + gradient.cwiseProduct(gradient);
    // denominator = adaDiag.cwiseSqrt().cwiseInverse();
    // beta = beta - eta * denominator.cwiseProduct(gradient);




    weight = weightiCpp(x, i, beta);
    logL[i] = logLikelihood(beta, y(i), weight);
    // Store the log-likelihood
    std::cout << logL[i] << std::endl;

    // std::cout << "y " << y(i) << std::endl;
    std::cout << "w " << weight << std::endl;
    // std::cout << "xrows " << x.row(i).rows() << " xcols " << x.row(i).cols() << std::endl;
    // std::cout << "brows " << beta.rows() << " bcols " << beta.cols() << std::endl;

    // Compute gradient on full data set

    // FIXME this is STILL not sparse.
    // gradient = gradSMCpp(zeros, x, i, beta, y(i), weight);
    // zeros = reset(zeros, x, i);
    // std::cout << "grows " << gradient.rows() << " gcols " << gradient.cols() << std::endl;
    // // Compute the squares of the elements of the gradient and use that for the hessian-ish-thing
    // adaDiag = adaDiag + gradient.cwiseProduct(gradient);
    // // std::cout << "arows " << adaDiag.rows() << " acols " << adaDiag.cols() << std::endl;
    //
    // // Add numerical tolerance, take 1/sqrt(A) to multiply
    // denominator = adaDiag.cwiseSqrt().cwiseInverse();
    // // std::cout << "drows " << denominator.rows() << " dcols " << denominator.cols() << std::endl;
    //
    // // Do the update
    // beta = beta - eta * denominator.cwiseProduct(gradient);
    // // std::cout << "brows " << beta.rows() << " bcols " << beta.cols() << std::endl;


  }


  return List::create(
    Named("adaDiag") = adaDiag,
    Named("beta") = beta,
    Named("logL") = logL
  );
}