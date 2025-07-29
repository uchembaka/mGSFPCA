// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
std::vector<Eigen::MatrixXd> get_SiCell(const Eigen::MatrixXd& data) {
  // Extract unique IDs
  Eigen::VectorXd IDs = data.col(0);
  std::vector<double> uniqueIDs;
  std::unordered_set<double> idSet;

  for (int i = 0; i < IDs.size(); ++i) {
    if (idSet.find(IDs[i]) == idSet.end()) {
      idSet.insert(IDs[i]);
      uniqueIDs.push_back(IDs[i]);
    }
  }
  int n = uniqueIDs.size();

  std::vector<Eigen::MatrixXd> SiCell(n);

  for (int i = 0; i < n; ++i) {
    // Find indices for current ID
    std::vector<int> ind;
    for (int j = 0; j < data.rows(); ++j) {
      if (data(j,0) == uniqueIDs[i]) {
        ind.push_back(j);
      }
    }
    int mi = ind.size();

    Eigen::VectorXd cYi(mi);
    for (int j = 0; j < mi; ++j) {
      cYi[j] = data(ind[j],3);
    }

    SiCell[i] = cYi * cYi.transpose();
  }

  return SiCell;
}
