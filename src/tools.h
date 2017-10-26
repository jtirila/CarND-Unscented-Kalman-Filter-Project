#include <vector>
#include "Eigen/Dense"


namespace Tools {
  /**
  * A helper function to convert from polar to cartesian coordinates.
  */
  Eigen::VectorXd ConvertToCartesian(const Eigen::VectorXd &polar_coords);

  /**
  * A helper function to convert from cartesian to polar coordinates.
  */
  Eigen::VectorXd ConvertToPolar(const Eigen::VectorXd &cart_coords);

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

};

