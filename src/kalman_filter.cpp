//**********************************************
// kalman_filter.cpp
//**********************************************
#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() 
{}

KalmanFilter::~KalmanFilter()
{}


void KalmanFilter::Predict() 
{
	//Predict state
	x_ = F_ *  x_;

	//Predict process noise.
	P_ = F_ * P_ * F_.transpose() + Q_;
  	
}

void KalmanFilter::Update(const VectorXd &z) 
{
	VectorXd y = z - (H_ * x_);

	MatrixXd S = H_ * P_ * H_.transpose() + R_;

	MatrixXd K = P_ * H_.transpose() * S.inverse();

	x_ = x_ + K * y;

	long x_size = x_.size();

	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	P_ = (I - K * H_) * P_;
}


void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
	
	float rho = sqrt(pow(x_[0],2) + pow(x_[1],2));
	float theta = atan2(x_[1],x_[0]);
	float rhodot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;

	VectorXd hx(3);

	 if (fabs(rho) < 0.0001)
	 {
		 rhodot = 0.001;
	 }
	 else
	 {
		 rhodot = (x_(0)*x_(2) + x_(1)*x_(3))/rho;
	 }

	hx << rho, theta, rhodot;

	VectorXd y = z - hx;

	//Angle normalization from udacity forums. Found
	//the formula from udacity forums.
	y[1] = atan2(sin(y[1]), cos(y[1]));

	MatrixXd S = H_ * P_ * H_.transpose() + R_;

	MatrixXd K = P_ * H_.transpose() * S.inverse();

	x_ = x_ + K * y;

	long x_size = x_.size();

	MatrixXd I = MatrixXd::Identity(x_size, x_size);

	P_ = (I - K * H_) * P_;
  	
}
