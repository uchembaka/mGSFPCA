// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "NLLK_Helper.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export(rng=false)]]
Eigen::VectorXd LSQg(const Eigen::VectorXd& cvec,
                     const List& SiCell,
                     const Eigen::MatrixXd& data,
                     int p, int k,
                     const Eigen::MatrixXd& B,
                     const Eigen::VectorXd& estGrid) {
  // Get dimensions
  int ncoefs = cvec.size();
  int M = B.rows();

  // Extract unique IDs
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

  // Create index boundaries
  Eigen::VectorXi indB(p+2);
  indB(0) = 1;
  for (int i = 1; i <= p; i++) {
    indB(i) = indB(i-1) + k;
  }

  // Split cvec into components
  Eigen::VectorXd CU = cvec.head(k*p);
  Eigen::VectorXd CDCSig2 = cvec.segment(indB(p)-1, cvec.size() - indB(p)+1);
  Eigen::VectorXd CD = CDCSig2.head(CDCSig2.size()-1);
  double Csig2 = CDCSig2(CDCSig2.size()-1);

  // Compute matrices
  Eigen::MatrixXd C = Map<Eigen::MatrixXd>(CU.data(), k, p);
  Eigen::MatrixXd U = B * C;
  Eigen::VectorXd D_vec = CD.array().exp();
  Eigen::MatrixXd D = D_vec.asDiagonal();
  double sig2 = std::exp(Csig2);

  // Gram-Schmidt orthogonalization
  List gs_result = gram_schmidt_with_deriv(U);
  Eigen::MatrixXd W = as<MatrixXd>(gs_result["W"]);
  std::vector<MatrixXd> DWDU = as<std::vector<MatrixXd>>(gs_result["DW"]);
  Eigen::MatrixXd Shat = W * D * W.transpose() + sig2 * MatrixXd::Identity(M, M);

  // Initialize gradient
  Eigen::VectorXd g = VectorXd::Zero(ncoefs);
  Eigen::MatrixXd DHT = D * W.transpose();
  Eigen::MatrixXd fullRes = MatrixXd::Zero(M, M);

  // Precompute ti and residuals for each unique_ID
  std::unordered_map<int, std::vector<int>> id_to_rows;
  for (int r = 0; r < data.rows(); ++r) {
    int id = data(r, 0);
    id_to_rows[id].push_back(r);
  }

  std::vector<VectorXi> ti_list(n);
  for (int i = 0; i < n; ++i) {
    const auto& pos = id_to_rows[unique_IDs(i)];
    Eigen::VectorXi ti(pos.size());
    for (size_t j = 0; j < pos.size(); ++j) {
      ti(j) = data(pos[j], 4) - 1; // Adjust index
    }
    // ti_list[i] = ti;

    // Compute residuals
    // const Eigen::VectorXi& ti = ti_list[i];
    Eigen::MatrixXd Sihat = Shat(ti, ti);
    // const Eigen::MatrixXd& Si = SiCell[i];
    Eigen::MatrixXd Si = as<MatrixXd>(SiCell[i]);
    Eigen::MatrixXd Res = Si - Sihat;
    fullRes(ti, ti) += Res;
  }

  // Compute gradient
  for (int j = 0; j < p; ++j) {
    Eigen::VectorXd gii = VectorXd::Zero(k);
    Eigen::MatrixXd dHdcl(M, p);
    dHdcl.setZero();

    for (int l = 0; l < k; ++l) {
      // Compute dHdcl for q >= j
      for (int q = j; q < p; ++q) {
        int idx = q * p + j;
        dHdcl.col(q) = DWDU[idx] * B.col(l);
      }
      Eigen::MatrixXd dHdclDHT = dHdcl * DHT;
      Eigen::MatrixXd dSihatdcl = dHdclDHT + dHdclDHT.transpose();
      gii(l) = 2.0 * (fullRes.cwiseProduct(-dSihatdcl)).sum();
    }
    g.segment(indB(j)-1, k) = gii;

    // Compute DDC and DDDc terms
    Eigen::VectorXd DDC = VectorXd::Zero(p);
    DDC(j) = D(j, j);
    Eigen::MatrixXd DDDc = -(W * DDC.asDiagonal() * W.transpose());
    g(indB(p) + j - 1) = 2.0 * (fullRes.cwiseProduct(DDDc)).sum();
  }

  // Compute sig2 term
  Eigen::MatrixXd sig2_eye = -sig2 * MatrixXd::Identity(M, M);
  g(ncoefs - 1) = 2.0 * (fullRes.cwiseProduct(sig2_eye)).sum();

  g /= (n * M);

  return g;
}
