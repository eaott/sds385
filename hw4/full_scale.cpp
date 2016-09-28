#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::Matrix2d;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::SparseVector;
using Eigen::RowMajor;
using namespace Rcpp;
typedef SparseVector<double> SVd;
typedef SparseMatrix<double, RowMajor> SMDR;
typedef Eigen::MappedSparseMatrix<double> MyMatrix;

/*
 * This file is almost identical to the other one.
 * Differences include setting m=1 everywhere for the LR
 * in order to not use an extra array. It also handles the sparsity.
 */

double TOL = 1e-6;
double EPSILON = 1e-8;

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
List test(MyMatrix xCWISE, VectorXd yVector) {
  /*
  * In this case, we need additional parameters to keep the running
  * values of beta and adaDiag going.
  *
  * Just go through the available values of x?
  */
  // SMDR x(xCWISE);
  // std::cout << "finished restructuring" << std::endl;
  std::cout << "got to c++" << std::endl;
  double eta = 1.0;
  double lambda = 0.01;



  const int N = xCWISE.cols(); // xCWISE is a P x N matrix!!!
  const int P = xCWISE.rows(); // xCWISE is a P x N matrix!!!
  VectorXd beta = VectorXd::Zero(P, 1);
  // VectorXd gradient = VectorXd::Zero(P, 1);
  // MyMatrix gradient(1, P); // column-oriented???
  SVd obs(P);
  VectorXd denominatorVector = VectorXd::Zero(P, 1);
  VectorXd adaDiagVec = VectorXd::Constant(P, 1, EPSILON);
  VectorXd gradientVec = VectorXd::Zero(P, 1);

  // TODO: this should work for the skip part.
  // I want skip == 1 to indicate not to do the update.
  // so when k == 0, skip == k - (-1) == 0 - (-1) == 1.
  VectorXi skipVec = VectorXi::Constant(P, 1, -1);


  double sumll = 0;
  for (int k=0; k < N; ++k) {
    // std::cout << k << std::endl;
    obs = xCWISE.innerVector(k);
    double y = yVector(k);
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
    double sum = obs.dot(beta);
    double expPart = exp(-1.0 * sum);
    double denominator = 1.0 + expPart;
    double weight = 1.0 / denominator;
    // update weight to be in [TOL, 1 - TOL].
    weight = (weight < TOL ? TOL : weight > 1 - TOL ? 1 - TOL : weight);


    ///////////////////////////////////////////
    // this section computes the log-likelihood
    ///////////////////////////////////////////
    double logL = y * log(weight) + (1 - y) * log(1 - weight);
    ///////////////////////////////////////////
    // this section computes the gradient
    ///////////////////////////////////////////
    double weightPart = y - weight;

    for (SVd::InnerIterator it(obs); it; ++it) {
      int index = it.index();

      // If I'm very very lucky, this should handle the
      // updates to the past L2 regularization iterations
      // that weren't actually performed.
      int skip = k - skipVec[index];
      double G_ii = adaDiagVec[index];
      double b_i = beta[index];
      // TODO without this for loop, takes ~11 seconds.
      // with it, it takes a long time. Need an approximation here...
      // for (int t = 0; t < skip - 1; t++) {
      //   G_ii += 4 * lambda * b_i * b_i;
      //   b_i = b_i * (1 - 2 * eta * lambda / sqrt(G_ii));
      // }
      adaDiagVec[index] = G_ii;
      beta[index] = b_i;

      double temp = -weightPart * it.value() + 2 * lambda * b_i;
      gradientVec[index] = temp;
      adaDiagVec[index] += temp * temp;
      double denom = 1.0 / sqrt(adaDiagVec[index]);

      beta[index] += -eta * denom * temp;


      // so here's the deal with the L2.
      // If we don't have it, then beta is already correct.
      // If we do, then each iteration, that element of beta
      // should be updated each iteration. We can prevent that
      // because we know what /would/ have been added to
      // adaDiagVec[i] based on the last version of gradientVec[i].
      // if skip == 1, then don't bother with anything.
      // See "Ridge SGD Sparse Weight Update" on
      // https://bryantravissmith.com/category/machine-learning/supervised-learning/classification/logistic-regression/
      // I think this looks like the following:
      // G_ii' = G_ii + gradF_i^2
      //       = G_ii + 4 * lambda * beta_i^2
      // beta_i' = beta_i * (1 - 2 * lambda * eta / sqrt(G_ii'))



      skipVec[index] = k;
    }
  }
  return List::create(
        Named("sumll") = sumll
      );
}









// // [[Rcpp::export]]
// List adagrad_sparse_sgd(SparseMatrix<double> xCWISE,
//                         VectorXd y,
//                         VectorXd beta,
//                         VectorXd adaDiag,
//                         double eta) {
//   /*
//   * In this case, we need additional parameters to keep the running
//   * values of beta and adaDiag going.
//   *
//   * Just go through the available values of x?
//   */
//   SMDR x(xCWISE);
//   std::cout << "finished restructuring" << std::endl;
//   const int N = x.rows();
//   const int P = x.cols();
//   if (y.rows() != N || beta.rows() != P || adaDiag.rows() != P) {
//     std::cout << "Wrong dimensions" << std::endl;
//   }
//   // Create matrices here.
//   VectorXd zeros = VectorXd::Zero(P, 1);
//   VectorXd gradient = VectorXd::Zero(P, 1);
//   VectorXd denominator = VectorXd::Zero(P, 1);
//   double weight = 0;
//
//   VectorXd logL = VectorXd::Zero(N, 1);
//   //
//   //   for (int k=0; k<x.outerSize(); ++k) {
//   //     for (SparseMatrix<double>::InnerIterator it(x,k); it; ++it)
//   //     {
//   //       it.value();
//   //       it.row();   // row index
//   //       it.col();   // col index (here it is equal to k)
//   //       it.index(); // inner index, here it is equal to it.row()
//   //     }
//   //   }
//
//
//   for (int i = 0; i < 1; i++) {
//
//
//     // Old way that kinda works if I use MatrixXd x.
//     // VectorXd row = x.row(i);
//     // weight = weightVCpp(row, beta);
//     // gradient = gradCpp(row, beta, y(i), weight);
//     // adaDiag = adaDiag + gradient.cwiseProduct(gradient);
//     // denominator = adaDiag.cwiseSqrt().cwiseInverse();
//     // beta = beta - eta * denominator.cwiseProduct(gradient);
//
//
//
//
//     weight = weightiCpp(x, i, beta);
//     logL[i] = logLikelihood(beta, y(i), weight);
//     // Store the log-likelihood
//     std::cout << logL[i] << std::endl;
//
//     // std::cout << "y " << y(i) << std::endl;
//     std::cout << "w " << weight << std::endl;
//     // std::cout << "xrows " << x.row(i).rows() << " xcols " << x.row(i).cols() << std::endl;
//     // std::cout << "brows " << beta.rows() << " bcols " << beta.cols() << std::endl;
//
//     // Compute gradient on full data set
//
//     // FIXME this is STILL not sparse.
//     // gradient = gradSMCpp(zeros, x, i, beta, y(i), weight);
//     // zeros = reset(zeros, x, i);
//     // std::cout << "grows " << gradient.rows() << " gcols " << gradient.cols() << std::endl;
//     // // Compute the squares of the elements of the gradient and use that for the hessian-ish-thing
//     // adaDiag = adaDiag + gradient.cwiseProduct(gradient);
//     // // std::cout << "arows " << adaDiag.rows() << " acols " << adaDiag.cols() << std::endl;
//     //
//     // // Add numerical tolerance, take 1/sqrt(A) to multiply
//     // denominator = adaDiag.cwiseSqrt().cwiseInverse();
//     // // std::cout << "drows " << denominator.rows() << " dcols " << denominator.cols() << std::endl;
//     //
//     // // Do the update
//     // beta = beta - eta * denominator.cwiseProduct(gradient);
//     // // std::cout << "brows " << beta.rows() << " bcols " << beta.cols() << std::endl;
//
//
//   }
//
//
//   return List::create(
//     Named("adaDiag") = adaDiag,
//     Named("beta") = beta,
//     Named("logL") = logL
//   );
// }