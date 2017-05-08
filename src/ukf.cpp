#include "ukf.h"
#include "tools.h"
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
  std_a_ = 8; // changed from 30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1; // changed from 30

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // Just constructed and hence not yet initialised
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension includes process noise
  n_aug_ = 7;

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //time initialise to zero
  time_us_ = 0.0;

  //initialise weights to be used for sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
  *  Initialization
  ****************************************************************************/

	if (!is_initialized_) {

		// Initialize state with first measurement

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

			double rho = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			double rhodot = meas_package.raw_measurements_(2);

			double px = rho * cos(phi);
			double py = rho * sin(phi);

			x_ << px, py, rho, phi, rhodot;

		}
		else {

			x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
		}
		// Initialize Covariance matrix 
		P_ << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
			-0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
			0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
			-0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
			-0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		time_us_ = meas_package.timestamp_;
		return;
	}

	else {
		double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
		Prediction(delta_t);
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			UpdateRadar(meas_package);
		}
		else {
			UpdateLidar(meas_package);
		}
		time_us_ = meas_package.timestamp_;
	}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
   
  /*****************************************************************************
  *  Prediction - Part 1 - Generate Sigma Points for k (current time)
  ****************************************************************************/
   //create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);
	x_aug.head(n_x_) = x_;
	x_aug.tail(n_aug_ - n_x_).setZero();

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
	P_aug.setZero();
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	//add noise
	P_aug.bottomRightCorner(n_aug_-n_x_, n_aug_ - n_x_) <<  std_a_*std_a_, 0,
															0, std_yawdd_*std_yawdd_;
	//create square root matrix
	MatrixXd A = P_aug.llt().matrixL();

	//create augmented sigma points
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
	//set first column of sigma point matrix
	Xsig_aug.col(0) = x_aug;
	//set remaining sigma points
	for (int i = 0; i < n_aug_; i++)
	{
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
	}

	/*****************************************************************************
	*  Prediction - Part 2 - Predict Sigma Points with noise for k+1 th time
	****************************************************************************/

	//predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x		= Xsig_aug(0, i);
		double p_y		= Xsig_aug(1, i);
		double v		= Xsig_aug(2, i);
		double yaw		= Xsig_aug(3, i);
		double yawd		= Xsig_aug(4, i);
		double nu_a     = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;

		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;

	}

	/*****************************************************************************
	*  Prediction - Part 3 - Predict Mean state and state Covariance of k+1 th time
	****************************************************************************/

	// Prediction of Mean And Covariance Matrix using predicted sigma points above
	double weight_0 = lambda_ / (lambda_ + n_aug_);
	weights_(0) = weight_0;
	for (int i = 1; i<2 * n_aug_ + 1; i++) {  //2n+1 weights
		double weight = 0.5 / (lambda_ + n_aug_);
		weights_(i) = weight;
	}

	//predict state mean
	x_.fill(0.0);
	for (size_t i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_ + weights_(i) * Xsig_pred_.col(i);
	}
	//predict state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

											   // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalizationz
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
		//std::cout << "P_" << i  << "=" << P << std::endl;
	}

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

  /*****************************************************************************
  *  UPDATE - Part 1 - Predict Measurement
  ****************************************************************************/

  //set measurement dimension, lidar can measure only px and py
	int n_z = 2;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
												// extract values for better readibility
		// measurement model
		Zsig(0, i) = Xsig_pred_(0, i);           //px
		Zsig(1, i) = Xsig_pred_(1, i);           //py
	}

	//mean predicted measurement
	//VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	//MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
												//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;

	S = S + R;

	/*****************************************************************************
	*  UPDATE - Part 2 - Update state
	****************************************************************************/

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (size_t i = 0; i < 2 * n_aug_ + 1; i++) {
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}
	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//create example vector for incoming radar measurement
	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_(0),
		meas_package.raw_measurements_(1);

	//residual
	VectorXd z_diff = z - z_pred;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	/*****************************************************************************
	*  UPDATE - Part 3 - Calculate Normalised Innovation Squared (NIS)
	****************************************************************************/
	NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;

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

  /*****************************************************************************
  *  UPDATE - Part 1 - Predict Measurement
  ****************************************************************************/

    //set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
	    // extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v   = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		Zsig(1, i) = atan2(p_y, p_x);                                 //phi
		Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//mean predicted measurement
	//VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		z_pred = z_pred + weights_(i) * Zsig.col(i);
	}

	//measurement covariance matrix S
	//MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											   //residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z, n_z);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	S = S + R;

  /*****************************************************************************
  *  UPDATE - Part 2 - Update state
  ****************************************************************************/

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);

	//calculate cross correlation matrix
	Tc.fill(0.0);
	for (size_t i = 0; i < 2 * n_aug_ + 1; i++) {
		//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}
	//calculate Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	//create example vector for incoming radar measurement
	VectorXd z = VectorXd(n_z);
	z << meas_package.raw_measurements_(0),
		meas_package.raw_measurements_(1),
		meas_package.raw_measurements_(2);

	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

	//update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K*S*K.transpose();

	/*****************************************************************************
	*  UPDATE - Part 3 - Calculate Normalised Innovation Squared (NIS)
	****************************************************************************/
	NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;
}

