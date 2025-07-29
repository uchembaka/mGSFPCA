// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "NLLK_Helper.h"

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export(rng=false)]]
double LSQf(const Eigen::VectorXd& cvec,
            const List& SiCell,
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
  Eigen::VectorXi unique_IDs = Eigen::Map<Eigen::VectorXi>(unique_IDs_vec.data(), unique_IDs_vec.size());
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

  Eigen::MatrixXd C = Map<MatrixXd>(CU.data(), k, p);
  Eigen::MatrixXd U = B * C;
  Eigen::VectorXd D_vec = CD.array().exp();
  Eigen::MatrixXd D = D_vec.asDiagonal();
  double sig2 = std::exp(Csig2);

  List gs_result = gram_schmidt_with_deriv(U);
  Eigen::MatrixXd W = as<MatrixXd>(gs_result["W"]);
  Eigen::MatrixXd Shat = W * D * W.transpose() + sig2 * MatrixXd::Identity(M, M);

  double f = 0.0;
  Eigen::VectorXd g = VectorXd::Zero(ncoefs);
  Eigen::MatrixXd DHT = D * W.transpose();
  Eigen::MatrixXd HD = W * D;


  for (int i = 0; i < n; i++) {
    std::vector<int> pos;
    for (int r = 0; r < data.rows(); r++) {
      if (static_cast<int>(data(r, 0)) == unique_IDs(i)) {
        pos.push_back(r);
      }
    }
    Eigen::VectorXi ti(pos.size());
    for (size_t j = 0; j < pos.size(); j++) {
      ti(j) = static_cast<int>(data(pos[j], 4)) - 1;
    }
    int mi = pos.size();

    Eigen::MatrixXd Sihat = Shat(ti, ti);
    Eigen::MatrixXd Si = as<MatrixXd>(SiCell[i]);
    Eigen::MatrixXd Res = Si - Sihat;
    double fi = Res.squaredNorm() / mi;
    f += fi;
  }
  f /= n;
  return f;
}
