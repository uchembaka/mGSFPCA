#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// Positive definite matrix correction
// [[Rcpp::export]]
Eigen::MatrixXd topdm_cpp(const Eigen::MatrixXd& sig) {
  const double EPS = 1e-6;
  const double ZERO = 1e-10;
  Eigen::MatrixXd sigma = sig;

  Eigen::LLT<Eigen::MatrixXd> llt(sigma);
  if (llt.info() == Eigen::NumericalIssue) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(sigma, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd d = svd.singularValues();

    for (int i = 0; i < d.size(); ++i) {
      if (d(i) <= ZERO) d(i) = EPS;
    }

    sigma = svd.matrixU() * d.asDiagonal() * svd.matrixV().transpose();
    sigma = 0.5 * (sigma + sigma.transpose());
  }
  return sigma;
}
