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
List test(MyMatrix xCWISE, VectorXd yVector, VectorXd training, const int nIter, const double lambda) {
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

  // Split the data into training and test data.
  SVd sparseTraining = training.sparseView();
  int nTrain = sparseTraining.sum();

  const int N = xCWISE.cols(); // xCWISE is a P x N matrix!!!
  SVd sparseTesting = (VectorXd::Constant(N, 1, 1) - training).sparseView();


  const int P = xCWISE.rows(); // xCWISE is a P x N matrix!!!
  VectorXd beta = VectorXd::Zero(P, 1);
  SVd obs(P);
  VectorXd adaDiagVec = VectorXd::Constant(P, 1, EPSILON);
  VectorXd logLVec = VectorXd::Zero(nIter * nTrain, 1);
  VectorXd logLTrainingVec = VectorXd::Zero(nIter * (N - nTrain), 1);
  VectorXi skipVec = VectorXi::Constant(P, 1, -1);


  double sumll = 0;
  for (int iter = 0; iter < nIter; ++iter) {
    int innerIndex = 0;
    const int iterCount = iter * nTrain;

    // Compute logL for the testing set.
    for (SVd::InnerIterator tr(sparseTesting); tr; ++tr) {
      int k = tr.index();
      // std::cout << k << std::endl;
      obs = xCWISE.innerVector(k);
      const double y = yVector(k);
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
      const double sum = obs.dot(beta);
      const double expPart = exp(-sum);
      const double denominator = 1.0 + expPart;
      double weight = 1.0 / denominator;
      // update weight to be in [TOL, 1 - TOL].
      weight = (weight < TOL) ?
      TOL :
        ((weight > 1 - TOL) ?
           1 - TOL :
           weight);


      ///////////////////////////////////////////
      // this section computes the log-likelihood
      ///////////////////////////////////////////
      // FIXME log-likelihood needs the lambda * ||beta||^2 term
      // that may be best maintained by keeping a running value?
      // since it's \Sum_i beta[i]^2, we know how to updated it
      // each iteration of the inner-most loop. Whatever the old value
      // was, subtract that (squared) and add the new value (squared).
      logLTrainingVec[iterCount + innerIndex] = y * log(weight) + (1 - y) * log(1 - weight);


    }
    for (SVd::InnerIterator tr(sparseTraining); tr; ++tr) {
      int k = tr.index();
      // std::cout << k << std::endl;
      obs = xCWISE.innerVector(k);
      const double y = yVector(k);
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
      const double sum = obs.dot(beta);
      const double expPart = exp(-sum);
      const double denominator = 1.0 + expPart;
      double weight = 1.0 / denominator;
      // update weight to be in [TOL, 1 - TOL].
      weight = (weight < TOL) ?
          TOL :
          ((weight > 1 - TOL) ?
             1 - TOL :
             weight);


      ///////////////////////////////////////////
      // this section computes the log-likelihood
      ///////////////////////////////////////////
      // FIXME log-likelihood needs the lambda * ||beta||^2 term
      // that may be best maintained by keeping a running value?
      // since it's \Sum_i beta[i]^2, we know how to updated it
      // each iteration of the inner-most loop. Whatever the old value
      // was, subtract that (squared) and add the new value (squared).
      logLVec[iterCount + innerIndex] = y * log(weight) + (1 - y) * log(1 - weight);


      ///////////////////////////////////////////
      // this section computes the gradient
      ///////////////////////////////////////////
      const double weightPart = y - weight;

      for (SVd::InnerIterator it(obs); it; ++it) {
        const int index = it.index();

        // If I'm very very lucky, this should handle the
        // updates to the past L2 regularization iterations
        // that weren't actually performed.
        const int skip = iterCount + innerIndex - skipVec[index];
        double G_ii = adaDiagVec[index];
        double b_i = beta[index];

        // Approx:
        // // Pretend like beta didn't change throughout the last `skip`
        // // iterations.
        // G_ii += (skip - 1) * 4 * lambda * lambda * b_i * b_i;
        // // Should this jump be bigger? According to the bryantravissmith
        // // website, maybe take fast_pow(1 - 2 * eta * lambda / sqrt(G_ii), (skip - 1))???
        // b_i = b_i * (1 - 2 * eta * lambda / sqrt(G_ii));
        G_ii = G_ii + (skip - 1) * 4 * lambda * lambda * b_i * b_i;
        b_i = b_i * (1 - 2 * eta * lambda / sqrt(G_ii));

        const double grad = -weightPart * it.value() + 2 * lambda * b_i;
        const double adaDiag = G_ii + grad * grad;
        // gradientVec[index] = temp;
        const double denom = 1.0 / sqrt(adaDiag);

        beta[index] = b_i - eta * denom * grad;
        adaDiagVec[index] = adaDiag;

        skipVec[index] = iterCount + innerIndex;
      }
      innerIndex++;
    }
  }
  const int curIndex = nIter * nTrain - 1;
  for (int j = 0; j < P; j++) {
    const int skip = curIndex - skipVec[j];
    double G_ii = adaDiagVec[j];
    double b_i = beta[j];
    G_ii = G_ii + (skip - 1) * 4 * lambda * lambda * b_i * b_i;
    if (skip >= 1) {
      b_i = b_i * (1 - 2 * eta * lambda / sqrt(G_ii));
    }
    beta[j] = b_i;
  }

  return List::create(
        Named("sumll") = sumll,
        Named("beta") = beta,
        Named("logL") = logLVec,
        Named("testLogL") = logLTrainingVec
      );
}
