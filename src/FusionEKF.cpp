#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  
  H_laser_ << 1,0,0,0,
             0,1,0,0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float d = measurement_pack.raw_measurements_[0];
      float th = measurement_pack.raw_measurements_[1];
      float v = measurement_pack.raw_measurements_[2];
      
      float py = d*std::sin(th);
      float px = d*std::cos(th);
      float vx = d*v/(py*std::tan(th)+px);
      float vy = d*v*tan(th)/(py*std::tan(th)+px);

      ekf_.x_ << px,py,vx,vy;
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER){
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
      ekf_.P_ = MatrixXd(4,4);
      ekf_.P_ << 1,0,0,0,
                 0,1,0,0,
                 0,0,1000,0,
                 0,0,0,1000;
    
      ekf_.F_ = MatrixXd(4,4);
      ekf_.F_ << 1,0,1,0,
                 0,1,0,1,
                 0,0,1,0,
                 0,0,0,1;
 
      previous_timestamp_ = measurement_pack.timestamp_;
      // done initializing, no need to predict or update
      is_initialized_ = true;
      return;
    }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  //calculate dt
  double dt = ((double)measurement_pack.timestamp_ - (double)previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  //update F
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
             
  //update Q
  double dt2 = dt * dt;
  double dt3 = dt * dt2;
  double dt4 = dt2 * dt2;
  double noise_ax = 9;
  double noise_ay = 9;
  
  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << dt4*noise_ax/4,      0       ,dt3*noise_ax/2,      0       ,
                   0       ,dt4*noise_ay/4,      0       ,dt3*noise_ay/2,
             dt3*noise_ax/2,      0       , dt2*noise_ax ,      0       ,
                   0       ,dt3*noise_ay/2,      0       , dt2*noise_ay ;
  
  ekf_.Predict();
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    
    //update R
    ekf_.R_ = R_radar_; 
    
    //update H
    Hj_ = tools.CalculateJacobian(ekf_.x_); 
    ekf_.H_ = Hj_; 
   
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
   
    //update R
    ekf_.R_ = R_laser_;
    
    //update H
    ekf_.H_ = H_laser_;
    
    ekf_.Update(measurement_pack.raw_measurements_);
  }

//    cout << "x_ = " << ekf_.x_ << endl;
//    cout << "P_ = " << ekf_.P_ << endl;
}
