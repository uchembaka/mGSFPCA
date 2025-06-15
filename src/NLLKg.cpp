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

  std::vector<int> unique_IDs_vec;
  std::unordered_set<int> seen;
  for (int r = 0; r < data.rows(); r++) {
    int id = static_cast<int>(data(r, 0));
    if (seen.find(id) == seen.end()) {
      seen.insert(id);
      unique_IDs_vec.push_back(id);
    }
  }
  Eigen::VectorXi unique_IDs = Eigen::Map<VectorXi>(unique_IDs_vec.data(), unique_IDs_vec.size());


  int n = unique_IDs.size();
  int M = B.rows();

  Eigen::VectorXi indB(p+2);
  indB(0) = 1;
  for (int i = 1; i <= p; i++) {
    indB(i) = indB(i-1) + k;
  }

  Eigen::VectorXd CU = cvec.head(k*p);
  Eigen::VectorXd CDCSig2 = cvec.segment(indB(p)-1, cvec.size() - indB(p)+1);
  Eigen::VectorXd CD = CDCSig2.head(CDCSig2.size()-1);
  double Csig2 = CDCSig2(CDCSig2.size()-1);

  Eigen::MatrixXd C = Map<Eigen::MatrixXd>(CU.data(), k, p);
  Eigen::MatrixXd U = B * C;
  Eigen::VectorXd D_vec = CD.array().exp();
  Eigen::MatrixXd D = D_vec.asDiagonal();
  double sig2 = std::exp(Csig2);

  List gs_result = gram_schmidt_with_deriv(U);
  Eigen::MatrixXd W = gs_result["W"];
  std::vector<MatrixXd> DWDU = gs_result["DW"];
  Eigen::MatrixXd Shat = W * D * W.transpose() + sig2 * MatrixXd::Identity(M, M);

  Eigen::VectorXd g = VectorXd::Zero(ncoefs);
  Eigen::MatrixXd DHT = D * W.transpose();
  Eigen::MatrixXd HD = W * D;

  std::unordered_map<int, std::vector<int>> id_to_rows;
  for (int r = 0; r < data.rows(); ++r) {
    int id = data(r, 0);
    id_to_rows[id].push_back(r);
  }

  // Precompute ti for each unique_ID
  std::vector<VectorXi> ti_list(n);
  for (int i = 0; i < n; ++i) {
    const auto& pos = id_to_rows[unique_IDs(i)];
    Eigen::VectorXi ti(pos.size());
    for (size_t j = 0; j < pos.size(); ++j) {
      ti(j) = data(pos[j], 4) - 1; // Adjust index
    }
    ti_list[i] = ti;
  }

    for (int i = 0; i < n; ++i) {

      const Eigen::VectorXi& ti = ti_list[i];
      Eigen::MatrixXd Sihat = Shat(ti, ti);
      const Eigen::MatrixXd& Si = SiCell[i];
      Eigen::MatrixXd SihatInv = Sihat.llt().solve(MatrixXd::Identity(ti.size(), ti.size()));
      Eigen::MatrixXd SihatInv_Si = SihatInv * Si;
      Eigen::MatrixXd dHdcl(M, p);
      Eigen::MatrixXd dSihatdcl(M, M);
      Eigen::MatrixXd dSihatdcl_i(ti.size(), ti.size());


      for (int j = 0; j < p; ++j) {
        for (int l = 0; l < k; ++l) {
          // Reset dHdcl
          dHdcl.setZero();

          // Compute dHdcl for q >= j
          for (int q = j; q < p; ++q) {
            int idx = q * p + j;
            dHdcl.col(q) = DWDU[idx] * B.col(l);
          }

          dSihatdcl.noalias() = dHdcl * DHT + HD * dHdcl.transpose();
          dSihatdcl_i = dSihatdcl(ti, ti);
          Eigen::MatrixXd tmp_mat = SihatInv * dSihatdcl_i;
          g(indB(j) + l - 1) += -(tmp_mat * SihatInv_Si).trace() + tmp_mat.trace();
        }

        // Compute DDC and DDDc terms
        Eigen::VectorXd DDC = VectorXd::Zero(p);
        DDC(j) = D(j, j);
        Eigen::MatrixXd DDDc = W * DDC.asDiagonal() * W.transpose();
        Eigen::MatrixXd DDDc_i = DDDc(ti, ti);
        Eigen::MatrixXd tmp_mat2 = SihatInv * DDDc_i;
        g(indB(p) + j - 1) += -(tmp_mat2 * SihatInv_Si).trace() + tmp_mat2.trace();
      }

      // Compute sig2 term
      int mi = ti.size();
      Eigen::MatrixXd sig2_eye = sig2 * MatrixXd::Identity(mi, mi);
      Eigen::MatrixXd tmp_mat3 = SihatInv * sig2_eye;
      g(ncoefs - 1) += -(tmp_mat3 * SihatInv_Si).trace() + tmp_mat3.trace();
    }

  g /= n;
  return g;
}
