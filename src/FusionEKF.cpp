#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF()
{
	is_initialized_ = false;

  	previous_timestamp_ = 0;

  	// initializing matrices
 	R_laser_ = MatrixXd(2, 2);
  	R_radar_ = MatrixXd(3, 3);
  	H_laser_ = MatrixXd(2, 4);
  	Hj_ = MatrixXd(3, 4);

  	ekf_.Q_ = MatrixXd(4,4);
	ekf_.F_ = MatrixXd(4,4);
	ekf_.P_ = MatrixXd(4,4);
	ekf_.x_ = VectorXd(4);


	H_laser_ << 1, 0, 0, 0,
		 		0, 1, 0, 0;

  	//measurement covariance matrix - laser
  	R_laser_ << 0.0225, 0,
  		        0, 		0.0225;

  	//measurement covariance matrix - radar
  	R_radar_ << 0.09,  0,      0,
  	      	    0,     0.0009, 0,
  	      	    0,     0,      0.09;

	ekf_.F_ <<  1, 0, 1, 0,
				0, 1, 0, 1,
				0, 0, 1, 0,
				0, 0, 0, 1;

	ekf_.P_ <<  1, 0, 0, 0,
				0 ,1, 0, 0,
				0, 0, 1000, 0,
				0, 0, 0, 1000;


	//Assign variance for prediction noise as 9
	//based on the question.
  	noise_ax = 9;
	noise_ay = 9;

}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
	if (!is_initialized_)
	{
		Initialize(measurement_pack);
	}
	else
	{
		Predict(measurement_pack);
		Update(measurement_pack);
	}
}


void FusionEKF::Initialize(const MeasurementPackage &measurement_pack)
{
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	{
		float rho = measurement_pack.raw_measurements_[0];
		float phi = measurement_pack.raw_measurements_[1];
		float rhodot = measurement_pack.raw_measurements_[2];

		//The input is of format rho, phi, rhodot
		//Convert into cartesian coordinate system.
		float x = measurement_pack.raw_measurements_[0] * sin(phi);
		float y = measurement_pack.raw_measurements_[0] * cos(phi);

		//Discard vx and vy as these are initial values and we do not have
		//much data to make sense of vx and vy
		ekf_.x_ << x, y , 0, 0;
		previous_timestamp_ = measurement_pack.timestamp_;

	}
	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
	{
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
		previous_timestamp_ = measurement_pack.timestamp_;
	}

	is_initialized_ = true;
}

void FusionEKF::Predict(const MeasurementPackage &measurement_pack)
{
	  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	  previous_timestamp_ = measurement_pack.timestamp_;

	  float dt_2 = dt * dt;
	  float dt_3 = dt_2 * dt;
	  float dt_4 = dt_3 * dt;

	  ekf_.F_(0,2) = dt;
	  ekf_.F_(1,3) = dt;

	  //To calculate process noise we need Q.
	  ekf_.Q_ << 	dt_4/4*noise_ax, 	0, 					dt_3/2*noise_ax, 	0,
			   	    0, 					dt_4/4*noise_ay, 	0, 					dt_3/2*noise_ay,
				    dt_3/2*noise_ax, 	0, 					dt_2*noise_ax, 		0,
					0, 					dt_3/2*noise_ay, 	0, 					dt_2*noise_ay;

	  ekf_.Predict();
}


void FusionEKF::Update(const MeasurementPackage &measurement_pack)
{
	  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
	  {

		  //First calculate Hj which is jacobian of
		  //state matrix
		  ekf_.H_= tools.CalculateJacobian(ekf_.x_);
		  ekf_.R_ = R_radar_;

		  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

	  }
	  else
	  {
		  ekf_.H_ = H_laser_;
		  ekf_.R_ = R_laser_;

		  ekf_.Update(measurement_pack.raw_measurements_);

	  }
}
