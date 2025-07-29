// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "NLLK_Helper.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export(rng=false)]]
Eigen::VectorXd NLLKg(const Eigen::VectorXd& cvec,
                      const std::vector<Eigen::MatrixXd>& SiCell,
                      const Eigen::MatrixXd& data,
                      int p, int k,
                      const Eigen::MatrixXd& B,
                      const Eigen::VectorXd& estGrid) {
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

  // --- build indB for parameter blocks ---
  Eigen::VectorXi indB(p+2);
  indB(0) = 1;
  for (int i = 1; i <= p; ++i) {
    indB(i) = indB(i-1) + k;
  }

  // --- unpack cvec into C, D, sigma2 ---
  Eigen::VectorXd CU      = cvec.head(k * p);
  Eigen::VectorXd tail    = cvec.segment(indB(p)-1, cvec.size() - indB(p)+1);
  Eigen::VectorXd CD      = tail.head(tail.size() - 1).array().exp();
  double           sig2   = std::exp(tail.tail<1>()(0));

  Eigen::MatrixXd C = Map<const MatrixXd>(CU.data(), k, p);
  Eigen::MatrixXd U = B * C;

  // --- orthogonalize ---
  List          gs     = gram_schmidt_with_deriv(U);
  Eigen::MatrixXd W    = gs["W"];
  std::vector<MatrixXd> DWDU = gs["DW"];  // length = p*p (per original idx = q*p + j)

  // --- build HD and DHT ---
  Eigen::MatrixXd D    = CD.asDiagonal();
  Eigen::MatrixXd HD   = W * D;
  Eigen::MatrixXd DHT  = D * W.transpose();
  int M = W.rows();

  // --- precompute full dShat / d cvec for every coefficient ---
  std::vector<MatrixXd> dShat_all(ncoefs);
  // 1) derivatives w.r.t. eigenfunction coefs
  for (int a = 0; a < k*p; ++a) {
    int j = a / k;       // column in C
    int l = a % k;       // row in C
    MatrixXd dHdcl = MatrixXd::Zero(M, p);
    for (int q = j; q < p; ++q) {
      int idx = q * p + j;
      dHdcl.col(q) = DWDU[idx] * B.col(l);
    }
    dShat_all[a] = dHdcl * DHT + HD * dHdcl.transpose();
  }
  // 2) derivatives w.r.t. eigenvalues
  for (int j = 0; j < p; ++j) {
    int a = k*p + j;
    VectorXd DDC = VectorXd::Zero(p);
    DDC(j) = CD(j);
    dShat_all[a] = W * DDC.asDiagonal() * W.transpose();
  }
  // 3) derivative w.r.t. sigma2
  {
    int a = ncoefs - 1;
    dShat_all[a] = sig2 * MatrixXd::Identity(M, M);
  }

  // --- map each unique ID to its row indices ---
  std::unordered_map<int, std::vector<int>> id_to_rows;
  id_to_rows.reserve(n);
  for (int r = 0; r < data.rows(); ++r) {
    id_to_rows[ static_cast<int>(data(r,0)) ].push_back(r);
  }

  std::vector<VectorXi> ti_list(n);
  for (int i = 0; i < n; ++i) {
    auto& pos = id_to_rows[ unique_IDs_vec[i] ];
    VectorXi ti(pos.size());
    for (size_t j = 0; j < pos.size(); ++j) {
      ti(j) = static_cast<int>(data(pos[j],4)) - 1;
    }
    ti_list[i] = ti;
  }

  // --- initialize gradient and accumulate over subjects ---
  VectorXd g = VectorXd::Zero(ncoefs);
  for (int i = 0; i < n; ++i) {
    const VectorXi& ti = ti_list[i];
    MatrixXd Sihat    = (W * D * W.transpose() + sig2 * MatrixXd::Identity(M,M))(ti, ti);
    MatrixXd Si       = SiCell[i];
    MatrixXd SihatInv = Sihat.llt().solve(MatrixXd::Identity(ti.size(), ti.size()));
    MatrixXd SihatInv_Si = SihatInv * Si;

    // loop over all coefficients
    for (int a = 0; a < ncoefs; ++a) {
      MatrixXd dS_i = dShat_all[a](ti, ti);// subset the full derivative matrix
      MatrixXd tmp = SihatInv * dS_i;
      double tr1   = (tmp * SihatInv_Si).trace();
      double tr2   = tmp.trace();

      g(a) += -(tr1) + tr2;
    }
  }

  g /= double(n);
  return g;
}
