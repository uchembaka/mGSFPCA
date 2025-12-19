// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Weighted inner product
inline double weighted_inner_product(const VectorXd& x,
                                     const VectorXd& y,
                                     const VectorXd& weights) {
  return (x.array() * weights.array() * y.array()).sum();
}


// Weighted norm
inline double weighted_norm(const VectorXd& x,
                            const VectorXd& weights) {
  return std::sqrt(weighted_inner_product(x, x, weights));
}


// [[Rcpp::export]]
Eigen::MatrixXd mwGS(const Eigen::MatrixXd& U,
                     const Eigen::VectorXd& weights) {
  const int D = U.rows();
  const int d = U.cols();

  MatrixXd W = U; // Start with copy of U
  VectorXd w(D);

  for (int i = 0; i < d; ++i) {
    // Normalize current vector
    double norm_i = weighted_norm(W.col(i), weights);
    if (norm_i > 1e-10) {
      W.col(i) /= norm_i;
    } else {
      W.col(i).setZero();
      continue;
    }

    // Orthogonalize subsequent vectors against current one
    for (int j = i + 1; j < d; ++j) {
      double coeff = weighted_inner_product(W.col(i), W.col(j), weights);
      W.col(j).noalias() -= coeff * W.col(i);
    }
  }

  return W;
}
