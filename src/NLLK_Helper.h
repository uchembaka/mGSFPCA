#ifndef NLLK_HELPER_H
#define NLLK_HELPER_H

#include <RcppEigen.h>
#include <cmath>

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

// Weighted Gram–Schmidt with derivatives
inline List gram_schmidt_with_deriv(const MatrixXd& U,
                                    const VectorXd& weights) {
  const int D = U.rows();
  const int d = U.cols();

  const MatrixXd I    = MatrixXd::Identity(D, D);
  const MatrixXd Wmat = weights.asDiagonal();

  std::vector<VectorXd> w(d, VectorXd(D));
  std::vector<VectorXd> q(d, VectorXd(D));
  std::vector<double>   norm_w(d), inv_norm3(d);
  std::vector<MatrixXd> grads(d * d, MatrixXd::Zero(D, D));
  std::vector<MatrixXd> grads_normed(d * d, MatrixXd::Zero(D, D));

  // initialize q with U columns
  for (int i = 0; i < d; ++i) q[i] = U.col(i);

  // ---- first vector ----
  w[0] = q[0];

  double norm2_0 = w[0].transpose() * Wmat * w[0]; // ||w_0||_w^2
  norm_w[0]      = std::sqrt(norm2_0);             // ||w_0||_w
  inv_norm3[0]   = 1.0 / (norm2_0 * norm_w[0]);    // 1 / ||w_0||_w^3

  grads[0] = I;  // dw_0 / dq_0 = I

  VectorXd Ww0 = Wmat * w[0];
  MatrixXd normMat0 =
    I / norm_w[0] - (w[0] * Ww0.transpose()) * inv_norm3[0];
  grads_normed[0] = normMat0 * grads[0];

  // temporary row-vectors for shear terms
  RowVectorXd qiG(D), wjG(D);

  // ---- remaining vectors ----
  for (int i = 1; i < d; ++i) {
    // weighted GS orthogonalization
    VectorXd cur = q[i];
    for (int j = 0; j < i; ++j) {
      double denom = w[j].transpose() * Wmat * w[j];
      double coeff = w[j].dot(Wmat * q[i]) / denom;
      cur.noalias() -= coeff * w[j];
    }
    w[i] = cur;

    double norm2_i = w[i].transpose() * Wmat * w[i];
    norm_w[i]      = std::sqrt(norm2_i);
    inv_norm3[i]   = 1.0 / (norm2_i * norm_w[i]);

    // normalization matrix: I/||w_i||_w - w_i (W w_i)^T / ||w_i||_w^3
    VectorXd Wwi = Wmat * w[i];
    MatrixXd normMat =
      I / norm_w[i] - (w[i] * Wwi.transpose()) * inv_norm3[i];

    // ---- gradients ----
    for (int k = 0; k <= i; ++k) {
      int idx = k * d + i;
      MatrixXd G = MatrixXd::Zero(D, D);

      if (k == i) {
        // derivative of w_i wrt q_i:
        //   dw_i/dq_i = I - sum_{j<i} P_j, P_j weighted proj
        int prev = (i - 1) * d + (i - 1);

        // last projection P_{i-1}
        double denom = w[i-1].transpose() * Wmat * w[i-1];
        VectorXd Ww_prev = Wmat * w[i-1];
        MatrixXd projection =
          (w[i-1] * Ww_prev.transpose()) / denom;

        G = grads[prev] - projection * grads[prev];
      } else {
        // shear terms (adapted to weighted inner product)
        for (int j = 0; j < i; ++j) {
          if (j < k) continue;

          double wj2 = w[j].dot(Wmat * w[j]);

          double wj_qi = w[j].dot(Wmat * q[i]);


          // weighted projections for grads
          qiG.noalias() = q[i].transpose() * Wmat * grads[k * d + j];
          wjG.noalias() = w[j].transpose() * Wmat * grads[k * d + j];

          MatrixXd term1 = w[j] * qiG;
          MatrixXd term2 = (2.0 * wj_qi / wj2) * (w[j] * wjG);
          MatrixXd term3 = wj_qi * grads[k * d + j];

          G.noalias() -= (term1 - term2 + term3) / wj2;
        }
      }

      grads[idx]       = G;
      grads_normed[idx] = normMat * G;
    }
  }

  // build orthonormal basis columns φ_i = w_i / ||w_i||_w
  MatrixXd W(D, d);
  for (int i = 0; i < d; ++i) {
    W.col(i) = w[i] / norm_w[i];
  }

  return List::create(
    _["W"]  = W,
    _["DW"] = grads_normed
  );
}

// convenience unweighted version
inline List gram_schmidt_with_deriv(const MatrixXd& U) {
  VectorXd weights = VectorXd::Ones(U.rows());
  return gram_schmidt_with_deriv(U, weights);
}

#endif
