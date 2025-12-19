// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "NLLK_Helper.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export(rng=false)]]
Eigen::VectorXd LSQg(const Eigen::VectorXd& cvec,
                     const std::vector<Eigen::MatrixXd>& SiCell,
                     const Eigen::MatrixXd& data,
                     int p, int k,
                     const Eigen::MatrixXd& B,
                     const Eigen::VectorXd& estGrid,
                     const VectorXd& weights) {
  int ncoefs = cvec.size();

  // --- extract unique IDs ---
  std::vector<int> unique_IDs_vec;
  unique_IDs_vec.reserve(data.rows());
  {
    std::unordered_set<int> seen;
    for (int r = 0; r < data.rows(); ++r) {
      int id = static_cast<int>(data(r, 0));
      if (seen.insert(id).second) {
        unique_IDs_vec.push_back(id);
      }
    }
  }
  int n = unique_IDs_vec.size();

  // --- indB for parameter blocks ---
  Eigen::VectorXi indB(p+2);
  indB(0) = 1;
  for (int i = 1; i <= p; ++i) {
    indB(i) = indB(i-1) + k;
  }

  // --- unpack parameters ---
  Eigen::VectorXd CU   = cvec.head(k * p);
  Eigen::VectorXd tail = cvec.segment(indB(p)-1, cvec.size() - indB(p)+1);

  Eigen::VectorXd CD_log = tail.head(tail.size() - 1);
  double          Csig2  = tail.tail<1>()(0);

  Eigen::VectorXd D_vec = CD_log.array().exp();
  double          sig2  = std::exp(Csig2);

  Eigen::MatrixXd C = Map<const MatrixXd>(CU.data(), k, p);
  Eigen::MatrixXd U = B * C;

  // --- Gram-Schmidt ---
  List gs = gram_schmidt_with_deriv(U, weights);
  Eigen::MatrixXd W = gs["W"];
  std::vector<MatrixXd> DWDU = gs["DW"];

  // --- D, HD, DHT ---
  Eigen::MatrixXd D   = D_vec.asDiagonal();
  Eigen::MatrixXd HD  = W * D;
  Eigen::MatrixXd DHT = D * W.transpose();
  int M               = W.rows();

  // --- precompute dShat_all (corrected indexing) ---
  std::vector<MatrixXd> dShat_all(ncoefs);

  // 1) derivatives w.r.t. C
  for (int a = 0; a < k*p; ++a) {
    int j = a / k;
    int l = a % k;

    MatrixXd dHdcl = MatrixXd::Zero(M, p);
    for (int q = j; q < p; ++q) {
      int idx = j * p + q;   // correct!!
      dHdcl.col(q) = DWDU[idx] * B.col(l);
    }
    dShat_all[a] = dHdcl * DHT + HD * dHdcl.transpose();
  }

  // 2) derivatives w.r.t. D
  for (int j = 0; j < p; ++j) {
    int a = k*p + j;
    VectorXd DDC = VectorXd::Zero(p);
    DDC(j) = D_vec(j);
    dShat_all[a] = W * DDC.asDiagonal() * W.transpose();
  }

  // 3) derivatives w.r.t. sigma^2
  {
    int a = ncoefs - 1;
    dShat_all[a] = sig2 * MatrixXd::Identity(M, M);
  }

  // --- subject intervals ---
  std::unordered_map<int, std::vector<int>> id_to_rows;
  id_to_rows.reserve(n);
  for (int r = 0; r < data.rows(); ++r) {
    id_to_rows[ static_cast<int>(data(r,0)) ].push_back(r);
  }

  std::vector<VectorXi> ti_list(n);
  std::vector<int> mi_list(n);
  for (int i = 0; i < n; ++i) {
    auto &pos = id_to_rows[ unique_IDs_vec[i] ];
    int mi    = pos.size();

    mi_list[i] = mi;
    VectorXi ti(mi);
    for (int j = 0; j < mi; ++j)
      ti(j) = static_cast<int>(data(pos[j],4)) - 1;

    ti_list[i] = ti;
  }

  // --- accumulate gradient ---
  VectorXd g = VectorXd::Zero(ncoefs);

  for (int i = 0; i < n; ++i) {
    const VectorXi &ti = ti_list[i];
    int mi             = mi_list[i];

    MatrixXd Shat = W * D * W.transpose() + sig2 * MatrixXd::Identity(M,M);
    MatrixXd Sihat = Shat(ti, ti);
    MatrixXd Si    = SiCell[i];
    MatrixXd Res   = Si - Sihat;

    for (int a = 0; a < ncoefs; ++a) {
      MatrixXd dS_i = dShat_all[a](ti, ti);

      // LSQ gradient:
      g(a) += -2.0 * (Res.cwiseProduct(dS_i)).sum() / double(mi);
    }
  }

  g /= double(n);
  return g;
}
