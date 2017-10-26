#include <iostream>
#include "tools.h"

/*
 * Convert a vector in polar coordinates into the cartesian coordinate system
 *
 * The input vector polar_coords is a three-element VectorXd containing rho (radius), phi (angle in radians) and
 * rhodot (range rate)
 *
 * Returns the vector with cartesian coordinates
 */
Eigen::VectorXd Tools::ConvertToCartesian(const Eigen::VectorXd &polar_coords){
  Eigen::VectorXd cartesian_coords(4);
  double r = polar_coords(0);
  double phi = polar_coords(1);

  cartesian_coords << std::cos(phi) * r, std::sin(phi) * r, 0, 0;

  return cartesian_coords;
}


Eigen::VectorXd Tools::CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                              const std::vector<Eigen::VectorXd> &ground_truth){

  Eigen::VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if(estimations.empty() || estimations.size() != ground_truth.size()){
    std::cout << "Error!";
    return rmse;
  }

  // Accumulate the residuals
  for(int i=0; i < estimations.size(); ++i){
    Eigen::VectorXd residuals = estimations[i] - ground_truth[i];
    Eigen::VectorXd sqr_res = residuals.array() * residuals.array();
    rmse += sqr_res;
  }

  rmse /= estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}
