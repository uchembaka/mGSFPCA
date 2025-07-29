#ifndef NLLK_HELPER_H
#define NLLK_HELPER_H

#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;
using namespace Eigen;

// Gram-Schmidt
inline List gram_schmidt_with_deriv(const MatrixXd& U) {
  const int D = U.rows();
  const int d = U.cols();

  std::vector<VectorXd> w(d, VectorXd(D));
  std::vector<VectorXd> q(d, VectorXd(D));
  std::vector<double> norm_w(d), inv_norm3(d);
  std::vector<MatrixXd> grads(d * d, MatrixXd::Zero(D, D));
  std::vector<MatrixXd> grads_normed(d * d, MatrixXd::Zero(D, D));

  const MatrixXd I = MatrixXd::Identity(D, D);

  for (int i = 0; i < d; ++i) q[i] = U.col(i);

  w[0] = q[0];
  norm_w[0] = w[0].norm();
  inv_norm3[0] = 1.0 / (norm_w[0] * norm_w[0] * norm_w[0]);
  grads[0] = I;
  grads_normed[0] = I / norm_w[0] - (w[0] * w[0].transpose()) * inv_norm3[0];

  VectorXd cur(D);
  RowVectorXd qiG(D), wjG(D);

  // Iterate
  for (int i = 1; i < d; ++i) {
    // GS orthogonalization
    cur = q[i];
    for (int j = 0; j < i; ++j) {
      double denom = w[j].squaredNorm();
      double coeff = w[j].dot(q[i]) / denom;
      cur.noalias() -= coeff * w[j];
    }
    w[i] = cur;
    norm_w[i] = w[i].norm();
    inv_norm3[i] = 1.0 / (norm_w[i] * norm_w[i] * norm_w[i]);

    // normalization matrix
    MatrixXd normMat = I / norm_w[i] - (w[i] * w[i].transpose()) * inv_norm3[i];

    // Gradients
    for (int k = 0; k <= i; ++k) {
      int idx = k * d + i;
      MatrixXd G = MatrixXd::Zero(D, D);

      if (k == i) {
        // w_i derivative wrt itself
        int prev = (i - 1) * d + (i - 1);
        double denom = w[i-1].squaredNorm();
        G = grads[prev] - (w[i-1] * w[i-1].transpose()) / denom;
      } else {
        // shear terms inline D_ji*grads[kd+j]
        for (int j = 0; j < i; ++j) {
          if (j < k) continue;
          // compute wj • qi and wjᵀwj
          double wj2 = w[j].squaredNorm();
          double wj_qi = w[j].dot(q[i]);
          qiG.noalias() = q[i].transpose() * grads[k * d + j];
          wjG.noalias() = w[j].transpose() * grads[k * d + j];
          G.noalias() -= (
            w[j] * qiG
          - (2.0 * wj_qi / wj2) * (w[j] * wjG)
            + wj_qi * grads[k * d + j]
          ) / wj2;
        }
      }
      grads[idx] = G;
      grads_normed[idx].noalias() = normMat * G;
    }
  }

  // Build W and DW
  MatrixXd W(D, d);
  for (int i = 0; i < d; ++i) W.col(i) = w[i] / norm_w[i];

  return List::create(
    _["W"] = W,
    _["DW"] = grads_normed
  );
}

#endif
