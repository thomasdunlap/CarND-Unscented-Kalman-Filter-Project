#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // Don't init until first measurement
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Set dimension of augment
  n_aug_ = 7;

  // Spread parameter
  lambda_ = 0;

  // Matrix of sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Weights
  weights_ = VectorXd(2 * n_aug_ + 1);

  // Noise
  R_radar = MatrixXd(3, 3);
  R_laser = MatrixXd(2, 2);

  // Time init
  time_us_ = 0;

  // NIS init
  NIS_radar_ = 0;
  NIS_laser_ = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rhodot = meas_package.raw_measurements_(2);

      // polar coordinates to cartesian:
      //x = rho * cos(angle)
      //y = rho * sin(angle)
      // x_[2] can be adjusted
      x_ << rho * cos(phi), rho * sin(phi), 4, rhodot * cos(phi), rhodot * sin(phi);

      // State covariance matrix
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 0, std_radphi_, 0,
            0, 0, 0, 0, std_radphi_;

      R_radar << std_radr_*std_radr_, 0, 0,
                 0, std_radphi_*std_radphi_, 0,
                 0, 0, std_radrd_*std_radrd_;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Init state; can adjust x_[2:]
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 4, 0.5, 0.0;

      // State covariance
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

      // noise
      R_laser << std_laspx_*std_laspx_, 0,
                 0, std_laspy_*std_laspy_;
    }

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
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
  // Spreading parameter
  lambda_ = 3 - n_x_;
  // Initiate sigma point matrix
  MatrixXd Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  // Square root of P
  MatrixXd A_ = P_.llt().matrixL();

  // Set sigma points as columns
  Xsig_.col(0) = x_;
  for(int i = 0; i < n_x_; i++) {
    Xsig_.col(i+1) = x_ + std::sqrt(lambda_ + n_x_) * A_.col(i);
    Xsig_.col(i+1 + n_x_) = x_ - std::sqrt(lambda_ + n_x_) * A_.col(i);
  }
  // Spreading parameter for augmented
  lambda_ = 3 - n_aug_;

  // Augmented state vector
  VectorXd x_aug_ = VectorXd(7);
  // Augmented covariance
  MatrixXd P_aug_ = MatrixXd(7, 7);
  // Augmented sigma points
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug_.head(5) = x_;
  x_aug_[5] = 0;
  x_aug_[6] = 0;

  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;

  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_.bottomRightCorner(2, 2) = Q;

  MatrixXd A_aug = P_aug_.llt().matrixL();

  // Sigma points
  Xsig_aug_.col(0) = x_aug_;
  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1) = x_aug_ + std::sqrt(lambda_ + n_aug_) * A_aug.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - std::sqrt(lambda_ + n_aug_) * A_aug.col(i);
  }


  //set vectors for each part added to x
  VectorXd vec1 = VectorXd(5);
  VectorXd vec2 = VectorXd(5);

  // Predict sigma points
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd calc_col = Xsig_aug_.col(i);
    double px = calc_col(0);
    double py = calc_col(1);
    double v = calc_col(2);
    double yaw = calc_col(3);
    double yawd = calc_col(4);
    double v_aug = calc_col(5);
    double v_yawdd = calc_col(6);

    //original
    VectorXd orig = calc_col.head(5);

    if(yawd > .001) {
      // If yaw dot is not zero
      vec1 << (v / yawd)*(sin(yaw + yawd * delta_t) - sin(yaw)),
              (v / yawd)*(-cos(yaw + yawd * delta_t) + cos(yaw)),
              0,
              yawd * delta_t,
              0;
    } else {
      // If yaw dot is zero - avoid division by zero
      vec1 << v * cos(yaw) * delta_t,
              v * sin(yaw) * delta_t,
              0,
              yawd * delta_t,
              0;
    }


    vec2 << .5 * delta_t * delta_t * cos(yaw) * v_aug,
            .5 * delta_t * delta_t * sin(yaw) * v_aug,
            delta_t * v_aug,
            .5 * delta_t * delta_t * v_yawdd,
            delta_t * v_yawdd;

    // Assign predicted sigma points
    Xsig_pred_.col(i) << orig + vec1 + vec2;
  }

  // Predicted state vector
  VectorXd x_pred = VectorXd(n_x_);

  // Predicted covariance Matrix
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);

  x_pred.fill(0.0);
  P_pred.fill(0.0);

  for(int i = 0; i < 2 * n_aug_ + 1; i++) {

    // Weights
    if (i == 0) {
      weights_(i) = lambda_ / (lambda_ + n_aug_);
    } else {
      weights_(i) = .5 / (lambda_ + n_aug_);
    }

    // Predict mean state
    x_pred += weights_(i) * Xsig_pred_.col(i);
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // Predict state covariance matrix
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;

    // Angle normalization
    if (x_diff(3) > M_PI) {
      x_diff(3) -= 2. * M_PI;
    } else if (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
    }
    P_pred += weights_(i) * x_diff * x_diff.transpose();
  }

  x_ = x_pred;
  P_ = P_pred;

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // Measurement dimension (px and py)
  int n_z = 2;

  // Sigma points matrix in plane
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // Predicted mean
  VectorXd z_pred = VectorXd(n_z);

  // Covariance measures
  MatrixXd S = MatrixXd(n_z, n_z);

  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // Sigma point transformation
    VectorXd state_vec = Xsig_pred_.col(i);
    double px = state_vec(0);
    double py = state_vec(1);

    Zsig.col(i) << px,
                   py;

    // Sum mean predicted measurements
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Sum covariance measurements
  for (int i = 0; i < 1 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add noise
  S += R_laser;

  // State vector for measurements
  VectorXd z = VectorXd(n_z);

  double meas_px = meas_package.raw_measurements_(0);
  double meas_py = meas_package.raw_measurements_(1);

  z << meas_px,
       meas_py;

  // Cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // Calculate cross correlation
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Angle normalization
    if (x_diff(3) > M_PI) {
      x_diff(3) -= 2. * M_PI;
    } else if (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
    }

    VectorXd z_diff = Zsig.col(i) - z_pred;

    Tc += weights_(i) * x_diff * z_diff.transpose();

  }

  VectorXd z_diff = z - z_pred;
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // Update state and covariance
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // Measurement dimensions (r, phi, r_dot)
  int n_z = 3;

  // Sigma points
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // Predicted mean state
  VectorXd z_pred = VectorXd(n_z);

  // Covariance
  MatrixXd S = MatrixXd(n_z,n_z);

  // Fill everything with zeros
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);
  double rho = 0;
  double phi = 0;
  double rho_d = 0;

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // Transform sigma points
    VectorXd state_vec = Xsig_pred_.col(i);
    double px = state_vec(0);
    double py = state_vec(1);
    double v = state_vec(2);
    double yaw = state_vec(3);
    double yaw_d = state_vec(4);

    rho = sqrt(px*px + py*py);
    phi = atan2(py, px);
    rho_d = (px * cos(yaw) * v + py * sin(yaw) * v) / rho;

    Zsig.col(i) << rho,
                   phi,
                   rho_d;

    // Calculate predicted measurement
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Covariance matrix S calculations
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (z_diff(1) > M_PI) {
      z_diff(1) -= 2. * M_PI;
    } else if (z_diff(1) < - M_PI) {
      z_diff(1) += 2. * M_PI;
    }
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Noise
  S += R_radar;

  // Measurement state vector
  VectorXd z = VectorXd(n_z);

  double meas_rho = meas_package.raw_measurements_(0);
  double meas_phi = meas_package.raw_measurements_(1);
  double meas_rhod = meas_package.raw_measurements_(2);

  z << meas_rho,
       meas_phi,
       meas_rhod;

  // Cross correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  // Assign cross correlation values
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    if (x_diff(3) > M_PI) {
      x_diff(3) -= 2. * M_PI;
    } else if (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
    }

    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (z_diff(1) > M_PI) {
      z_diff(1) -= 2. * M_PI;
    } else if (z_diff(1) < -M_PI) {
      z_diff(1) += 2. * M_PI;
    }

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z_diff = z - z_pred;

  if (z_diff(1) > M_PI) {
    z_diff(1) -= 2. * M_PI;
  } else if (z_diff(1) < -M_PI) {
    z_diff(1) += 2. * M_PI;
  }

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // Update state and covariance
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

}
