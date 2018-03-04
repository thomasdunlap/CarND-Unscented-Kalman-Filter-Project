#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // Sum of squared residuals
  for (int i=0; i < estimations.size(); ++i) {
    VectorXd residuals = estimations.size() - ground_truth[i];
    residuals = residuals.array() * residuals.array();
    rmse += residuals;
  }

  // Mean
  rmse /= estimations.size();

  rmse = sqrt(rmse.array());

  return rmse;
}
