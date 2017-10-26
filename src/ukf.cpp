#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false; // TODO

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  lambda_ = 3 - n_x_;
  spreading_coeff_ = std::sqrt(lambda_ + n_x_);

  Tc_.fill(0.0);

  weights_.fill(1 / (2 * (lambda_ + Xsig_pred_.cols())));
  weights_(0) = lambda_ / (lambda_ + Xsig_pred_.cols());
  sum_weights_ = ((Xsig_pred_.cols() - 1) * weights_(1) + weights_(0));

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.


  Hint: one or more values initialized above might be wildly off...
  */

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */


  // For this project, just use a guard against weird rows in the measurement data. If these did really occur,
  // something more sophisticated would be needed.
  assert(measurement_pack.sensor_type_ == MeasurementPackage::RADAR
         || measurement_pack.sensor_type_ == MeasurementPackage::LASER);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR &&
      (M_PI < measurement_pack.raw_measurements_(1) || -M_PI > measurement_pack.raw_measurements_(1))) {

    // Do nothing with radar measurements where the measured angle is not between -PI and PI. This is probably
    // not necessary as one or two spurious measurement would do little harm in the big picture. However, I guess
    // it is one feasible strategy to deal with such values to just ignore them.
    cout << "Skipping a weird measurement, radar data with values " << measurement_pack.raw_measurements_ << " \n";
    return;
  }

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      VectorXd cart_coord = Tools::ConvertToCartesian(measurement_pack.raw_measurements_);

      x_(0) = cart_coord(0);
      x_(1) = cart_coord(1);
      // No need to set the rest.
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    // Debug:
    cout << "Initialization done\n";

    return;
  }

  // Start the processing if initialization had been performed previously already

  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && !use_laser_) {
    cout << "Ignoring laser measurements to investigate radar processing performance.\n";
    return;
  }

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && !use_radar_) {
    cout << "Ignoring radar measurements to investigate radar processing performance.\n";
    return;
  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/


  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  // Modify the F and Q matrices so that the time is integrated; these modifications will be in effect when entering the
  // Prediction and Update steps.

  Predict(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Call the function that will perform extended KF on the inputs
    UpdateRadar(measurement_pack);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Call the function that will perform basic KF on the inputs
    UpdateLidar(measurement_pack);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Predict(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  Augment();
  PredictSigmaPoints(delta_t);
  PredictMeanCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */



// void UKF::UpdateRadar(MeasurementPackage meas_package) {
//
//   for(int col = 0; col < Xsig_pred_.cols(); col++){
//     VectorXd xdiff = Xsig_.col(col) - x_;
//     VectorXd zdiff = Xsig_pred_.col(col) - z_pred_;
//     Tc_ += weights_(col) * xdiff * zdiff.transpose();
//   }
//
//   // Calculate Kalman gain K;
//   MatrixXd K = Tc_ * S.inverse();
//
//   // Update state mean and covariance matrix
//   x_ += K * (z_ - z_pred);
//   P_ -= K * S * K.transpose();
//
//
//
//
//   /**
//   TODO:
//
//   Complete this function! Use radar data to update the belief about the object's
//   position. Modify the state vector, x_, and covariance, P_.
//
//   You'll also need to calculate the radar NIS.
//   */
// }

void UKF::Augment(){
  x_aug_.head(n_x_) = x_;
  P_aug_.block(0, 0, n_x_, n_x_) = P_;
  P_aug_.block(0, n_x_, n_x_, 2) = MatrixXd::Zero(n_x_, 2);
  P_aug_.block(n_x_, 0, 2, n_x_) = MatrixXd::Zero(2, n_x_);
  P_aug_(n_x_, n_x_) = std_a_ * std_a_;
  P_aug_(n_x_ - 1, n_x_) = 0;
  P_aug_(n_x_, n_x_ - 1) = 0;
  P_aug_(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;


  MatrixXd A = P_aug_.llt().matrixL();

  MatrixXd x_aug_repl = x_aug_.replicate(1, n_aug_);

  Xsig_aug_.col(0) = x_aug_;
  Xsig_aug_.block(0, 1, n_aug_, n_aug_) = x_aug_repl + spreading_coeff_ * A;
  Xsig_aug_.block(0, n_aug_ + 1, n_aug_, n_aug_) = x_aug_repl - spreading_coeff_ * A;
}


/*
 * TODO: Generate sigma points. Maybe not necessary?
 */
void UKF::PredictSigmaPoints(double delta_t) {

  // Augmented state vector:
  // [px, py, v, phi, phidot, va, vphidotdot]
  VectorXd x_aug;
  VectorXd acc_vec(5);
  acc_vec.fill(0.0);
  VectorXd noise_vec(5);

  double px, py, v, phi, phidot, va, phidotdot;

  for (int column = 0; column < Xsig_aug_.cols(); column++) {
    px = Xsig_aug_(0, column);
    py = Xsig_aug_(1, column);
    v = Xsig_aug_(2, column);
    phi = Xsig_aug_(3, column);
    phidot = Xsig_aug_(4, column);
    va = Xsig_aug_(5, column);
    phidotdot = Xsig_aug_(6, column);


    if (std::abs(phidot) < 0.00001) {
      acc_vec(0) = 0.0;
      acc_vec(1) = 0.0;
    } else {
      acc_vec(0) = v / phidot * (std::sin(phi + phidot * delta_t) - std::sin(phi));
      acc_vec(1) = v / phidot * (-std::cos(phi + phidot * delta_t) + std::cos(phi));
    }
    // acc_vec(2) no change
    acc_vec(3) = phidot * delta_t;
    // acc_vec(4) no change


    double half_delta_t_sqr = delta_t * delta_t / 2;

    noise_vec(0) = half_delta_t_sqr * std::cos(phi) * va;
    noise_vec(1) = half_delta_t_sqr * std::sin(phi) * va;
    noise_vec(2) = delta_t * va;
    noise_vec(3) = half_delta_t_sqr * phidotdot;
    noise_vec(4) = delta_t * phidotdot;

    Xsig_pred_.col(column) = (VectorXd(5) << px, py, v, phi, phidot).finished() + acc_vec + noise_vec;
  }
}

void UKF::PredictMeanCovariance() {
  x_.fill(0.0);

  for (int col = 0; col < Xsig_pred_.cols(); col++) {
    x_ += weights_(col) * Xsig_pred_.col(col) / sum_weights_;
  }

  Eigen::VectorXd x_tmp;

  for (int col = 0; col < Xsig_pred_.cols(); col++) {
    x_tmp = Xsig_pred_.col(col);
    P_ += weights_(col) * (x_tmp - x_) * (x_tmp - x_).transpose() / sum_weights_;
  }
}


