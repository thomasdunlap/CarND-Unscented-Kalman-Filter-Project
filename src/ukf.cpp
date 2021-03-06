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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.6;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55;

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

  // State length
  n_x_ = 5;

  // Augmented state length
  n_aug_ = n_x_ + 2;

  // Lambda values to help initiate weights
  lambda_ = 3 - n_aug_;

  // Sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Set weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill((double) 0.5 / (n_aug_ + lambda_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // Initialization
  if (!is_initialized_) {

    // State vector
    x_.fill(0.0);

    // Covariance matrix
    P_ << 1.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0, 0.0,
          0.0, 0.0, 0.0, 0.0, 1.0;

    // Initatiate radar measurements
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);
      // x
      x_(0) = rho * cos(phi);
      // y
      x_(1) = rho * sin(phi);
      // v
      x_(2) = rho_dot;

    }

    // If laser
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

      // Adjust for x,y values too close to zero
      if (fabs(x_(0)) < 0.001) {
        x_(0) = 0.001;
      }
      if (fabs(x_(1)) < 0.001) {
        x_(1) = 0.001;
      }
    }
    // Set as initialized
    is_initialized_ = true;
    // Update time
    time_us_ = meas_package.timestamp_;
    return;
  }

    // Change in time - store for future use
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
    time_us_ = meas_package.timestamp_;

    // Predict locatization
    Prediction(dt);

    // Update measured values
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
/*******************************************************************************
 * Augmentation
 ******************************************************************************/

  // Augmented state vector
  VectorXd x_aug = VectorXd(7);

  // Augmented covariance matrix
  MatrixXd P_aug = MatrixXd(7, 7);

  // Augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Assign values
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;


  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // Square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  MatrixXd B = sqrt(lambda_ + n_aug_) * A;

  // Augmented sigma points
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.leftCols(n_aug_ + 1).rightCols(n_aug_) = x_aug.replicate(1, B.cols()) + B;
  Xsig_aug.rightCols(n_aug_) = x_aug.replicate(1, B.cols()) - B;

/*******************************************************************************
 * Prediction
 ******************************************************************************/

  // Predict sigma points
  for (int i = 0; i< 2 * n_aug_ + 1; i++) {
    // Better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    double delta_t2 = delta_t * delta_t;
    double yawd_dt = yawd * delta_t;

    // Control for zero-like division
    if (fabs(yawd) > 0.001) {
        Xsig_pred_(0, i) = p_x + v/yawd * (sin(yaw + yawd_dt) - sin(yaw)) +
                          (0.5 * delta_t2 * cos(yawd) * nu_a);
        Xsig_pred_(1, i) = p_y + v/yawd * (-cos(yaw + yawd_dt) + cos(yaw)) +
                          (0.5 * delta_t2 * sin(yawd) * nu_a);
    }
    else {
        Xsig_pred_(0, i) = p_x + v*delta_t * cos(yaw) +
                          (0.5*delta_t2 * cos(yawd) * nu_a);
        Xsig_pred_(1, i) = p_y + v*delta_t * sin(yaw) +
                          (0.5*delta_t2 * sin(yawd) * nu_a);
    }

    Xsig_pred_(2, i) = v + nu_a*delta_t;
    Xsig_pred_(3, i) = yaw + yawd_dt + 0.5*nu_yawdd * delta_t2;
    Xsig_pred_(4, i) = yawd + nu_yawdd * delta_t;
  }
/*******************************************************************************
 * Predicted Mean Covariance
 ******************************************************************************/

  // Predict state mean
  MatrixXd x_sum = Xsig_pred_.array().rowwise() * weights_.transpose().array();
  x_ = x_sum.rowwise().sum();

  // Predict state covariance matrix
  MatrixXd X_abs = Xsig_pred_.array().colwise() - x_.array();

  // Angle normalization
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
    X_abs(3, i) = NormalizeAngle(X_abs(3, i));
  }

  X_abs_w_ = X_abs.array().rowwise() * weights_.transpose().array();

  P_ = X_abs_w_ * X_abs.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

/*******************************************************************************
 * Predicted Lidar Measurement
 ******************************************************************************/

  // Sigma points
  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);

  // Predicted state
  VectorXd z_pred = VectorXd(2);

  // Measured state
  VectorXd z = VectorXd(2);

  // Measurement covariance
  MatrixXd S = MatrixXd(2, 2);

  // Sigma point transformation
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
      Zsig(0, i) = Xsig_pred_(0, i);
      Zsig(1, i) = Xsig_pred_(1, i);
  }

  // Predicted state
  MatrixXd z_sum = Zsig.array().rowwise() * weights_.transpose().array();
  z_pred = z_sum.rowwise().sum();

  // Calculations for S
  MatrixXd Z_abs = Zsig.array().colwise() - z_pred.array();
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
    // Normalize angles
    Z_abs(1, i) = NormalizeAngle(Z_abs(1, i));
  }
  MatrixXd Z_abs_w = Z_abs.array().rowwise() * weights_.transpose().array();
  S = Z_abs_w * Z_abs.transpose();

   // Measurement noise
  S(0, 0) += std_laspx_ * std_laspx_;
  S(1, 1) += std_laspy_ * std_laspy_;

/*******************************************************************************
 * UKF Radar Update
 ******************************************************************************/

  MatrixXd Tc = X_abs_w_ * Z_abs.transpose();

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  z = meas_package.raw_measurements_;

  // Difference between predicted and measured
  VectorXd z_diff = z - z_pred;

  // Normalize angle
  z_diff(1) = NormalizeAngle(z_diff(1));

  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  You'll also need to calculate the radar NIS.
  */
/*******************************************************************************
 * Predicted Radar Measurement
 ******************************************************************************/



  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(3);
  VectorXd z = VectorXd(3);

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(3, 3);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double y = Xsig_pred_(3, i);

      Zsig(0, i) = sqrt((px * px) + (py * py));
      Zsig(1, i) = atan2(py, px);
      Zsig(2, i) = (((px * cos(y)) + (py * sin(y))) * v) / Zsig(0, i);
  }

  //calculate mean predicted measurement
  MatrixXd z_sum = Zsig.array().rowwise() * weights_.transpose().array();
  z_pred = z_sum.rowwise().sum();

  //calculate measurement covariance matrix S
  MatrixXd Z_abs = Zsig.array().colwise() - z_pred.array();

  // Normalize angles
  for (int i=0; i < 2 * n_aug_ + 1; i++) {
    Z_abs(1, i) = NormalizeAngle(Z_abs(1, i));
  }

  MatrixXd Z_abs_w = Z_abs.array().rowwise() * weights_.transpose().array();

  S = Z_abs_w * Z_abs.transpose();

  // Add measurement noise covariance matrix
  S(0, 0) += std_radr_ * std_radr_;
  S(1, 1) += std_radphi_ * std_radphi_;
  S(2, 2) += std_radrd_ * std_radrd_;


/*******************************************************************************
 * UKF Radar Update
 ******************************************************************************/

  MatrixXd Tc = X_abs_w_ * Z_abs.transpose();

  // Kalman gain K
  MatrixXd K = Tc * S.inverse();

  z = meas_package.raw_measurements_;
  // Difference between predicted and measured
  VectorXd z_diff = z - z_pred;

  // Normalize angles
  z_diff(1) = NormalizeAngle(z_diff(1));

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  // Update state and covariance
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}

/*
  Normalizes angle if outside measurement space.
*/
double UKF::NormalizeAngle(double angle) {

  while (angle > M_PI) angle -= 2.0 * M_PI;
  while (angle < -M_PI) angle += 2.0 * M_PI;

  return angle;
}
