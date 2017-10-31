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
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.355;

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


  n_aug_ = 7;

  lambda_ = 3 - n_aug_;
  spreading_coeff_ = std::sqrt(lambda_ + n_aug_);

  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.fill(0.0);

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  weights_ = Eigen::VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);
  // FIXME: these are potentially off, check lectures
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  sum_weights_ = ((n_aug_ - 1) * weights_(1) + weights_(0));

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

      Eigen::VectorXd cart_coord = Tools::ConvertToCartesian(measurement_pack.raw_measurements_);

      x_(0) = cart_coord(0);
      x_(1) = cart_coord(1);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = measurement_pack.raw_measurements_(0);
      x_(1) = measurement_pack.raw_measurements_(1);
    }


    P_ = P_.Identity(P_.rows(), P_.cols());

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

  Predict(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  int n_z;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    n_z = 3;
  else
    n_z = 2;

  Eigen::MatrixXd S(n_z, n_z);
  Eigen::VectorXd z_pred(n_z);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  S.fill(0.0);
  z_pred.fill(0.0);
  Tc.fill(0.0);
  Zsig.fill(0.0);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR){

    PredictRadar(Zsig, z_pred, S);

    // Call the function that will perform extended KF on the inputs

  } else  {
    PredictLidar(Zsig, z_pred, S);
  }
  Update(measurement_pack, Zsig, S, z_pred, Tc);

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


void UKF::Update(MeasurementPackage meas_package, Eigen::MatrixXd& Zsig, Eigen::MatrixXd& S,
                      Eigen::VectorXd& z_pred, Eigen::MatrixXd& Tc) {



  for (int col = 0; col < Xsig_pred_.cols(); col++) {
    VectorXd xdiff = Xsig_pred_.col(col) - x_;
    xdiff(3) = Tools::NormalizeAngle(xdiff(3));
    VectorXd zdiff = Zsig.col(col) - z_pred;
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
      zdiff(1) = Tools::NormalizeAngle(zdiff(1));
    Tc += weights_(col) * xdiff * zdiff.transpose();
  }


  // Calculate Kalman gain K;
  Eigen::MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;

  // Update state mean and covariance matrix
  VectorXd zdiff = meas_package.raw_measurements_ - z_pred;
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    zdiff(1) = Tools::NormalizeAngle(zdiff(1));
  x_ += K * zdiff;
  P_ -= K * S * K.transpose();

  Tools::Calculate_Nis(S_inv, zdiff);


  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

void UKF::Augment(){

  Eigen::VectorXd x_aug = Eigen::VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;

  Eigen::MatrixXd P_aug = Eigen::MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  P_aug.block(0, 0, n_x_, n_x_) = P_;

  P_aug.block(0, n_x_, n_x_, 2) = MatrixXd::Zero(n_x_, 2);
  P_aug.block(n_x_, 0, 2, n_x_) = MatrixXd::Zero(2, n_x_);
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_, n_x_ + 1) = 0;
  P_aug(n_x_ + 1, n_x_) = 0;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  MatrixXd A = P_aug.llt().matrixL();

  MatrixXd x_aug_repl = x_aug.replicate(1, n_aug_);

  Xsig_aug_.col(0) = x_aug;
  Xsig_aug_.block(0, 1, n_aug_, n_aug_) = x_aug_repl + spreading_coeff_ * A;
  Xsig_aug_.block(0, n_aug_ + 1, n_aug_, n_aug_) = x_aug_repl - spreading_coeff_ * A;
}


void UKF::PredictSigmaPoints(double delta_t) {

  // Augmented state vector:
  // [px, py, v, phi, phidot, va, vphidotdot]
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


    if (std::fabs(phidot) < 0.001) {
      acc_vec(0) = v * std::cos(phi) * delta_t;
      acc_vec(1) = v * std::sin(phi) * delta_t;
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

/**
 *
 * @param Zsig
 * @param z_pred
 * @param S
 */
void UKF::PredictRadar(Eigen::MatrixXd& Zsig, Eigen::VectorXd& z_pred, Eigen::MatrixXd& S) {


  double px, py, v, phi, phidot;
  double rho, phi_meas, phidot_meas;

  for(int col = 0; col < Xsig_pred_.cols(); col++ ){
    px = Xsig_pred_(0, col);
    py = Xsig_pred_(1, col);
    v = Xsig_pred_(2, col);
    phi = Xsig_pred_(3, col);

    rho = std::sqrt(px * px + py * py);
    phi_meas = std::atan2(py, px);

    if(std::abs(rho) < 0.001)
      phidot_meas = 0;
    else
      phidot_meas = (px * std::cos(phi) * v + py * std::sin(phi) * v) / rho;

    Zsig(0, col) = rho;
    Zsig(1, col) = phi_meas;
    Zsig(2, col) = phidot_meas;

    z_pred += weights_(col) * Zsig.col(col);
  }

  Eigen::VectorXd diff;
  for(int col = 0; col < Zsig.cols(); col++){
    diff = Zsig.col(col) - z_pred;

    diff(1) = Tools::NormalizeAngle(diff(1));
    S += weights_(col) * diff * diff.transpose();

  }

  MatrixXd noise(3, 3);
  noise.fill(0.0);
  noise(0, 0) = std_radr_ * std_radr_;
  noise(1, 1) = std_radphi_ * std_radphi_;
  noise(2, 2) = std_radrd_ * std_radrd_;

  S += noise;

}


/**
 *
 * @param Zsig
 * @param z_pred
 * @param S
 */
void UKF::PredictLidar(Eigen::MatrixXd& Zsig, Eigen::VectorXd& z_pred, Eigen::MatrixXd& S) {


  for (int col = 0; col < Xsig_pred_.cols(); col++) {
    Zsig(0, col) = Xsig_pred_(0, col);
    Zsig(1, col) = Xsig_pred_(1, col);

    z_pred += weights_(col) * Zsig.col(col);
  }

  Eigen::VectorXd diff;
  for (int col = 0; col < Zsig.cols(); col++) {
    diff = Zsig.col(col) - z_pred;
    S += weights_(col) * diff * diff.transpose();

  }

  MatrixXd noise(2, 2);
  noise.fill(0.0);
  noise(0, 0) = std_laspx_;
  noise(1, 1) = std_laspy_;

  S += noise;
}

void UKF::PredictMeanCovariance() {
  x_.fill(0.0);

  for (int col = 0; col < Xsig_pred_.cols(); col++) {
    x_ += weights_(col) * Xsig_pred_.col(col);
  }

  P_.fill(0.0);
  Eigen::VectorXd xdiff;

  for (int col = 0; col < Xsig_pred_.cols(); col++) {
    xdiff = Xsig_pred_.col(col) - x_;
    xdiff(3) = Tools::NormalizeAngle(xdiff(3));
    P_ += weights_(col) * xdiff * xdiff.transpose();
  }
}


