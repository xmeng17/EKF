#include "kalman_filter.h"
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  //predict mean and variance
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

MatrixXd KalmanFilter::Gain() {
  //calculate Karman Gain
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  return K;
}

void KalmanFilter::Update(const VectorXd &z) {
  
  //calculate Karman Gain
  MatrixXd K = Gain();
 
  //update mean
  VectorXd y = z - H_ * x_;
  x_ = x_ + K * y;
  
  //update variance
  MatrixXd I = MatrixXd::Identity(x_.size(),x_.size());
  P_ = (I - K * H_) * P_; 
}

VectorXd KalmanFilter::StateToMeasurement(){
  
  VectorXd h(3);

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float d = std::sqrt(px * px + py * py);
  h << d,atan2(py,px),(px*vx+py*vy)/d;
  return h;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  //calculate Karman Gain
  MatrixXd K = Gain();
   
  //update mean
  VectorXd h = StateToMeasurement();
  VectorXd y = z - h;
  
  //normalize angle
  while (y(1)>M_PI) y(1) -= 2 * M_PI;
  while (y(1)<-M_PI) y(1) += 2 * M_PI;

  x_ = x_ + K * y;
  
  //update variance
  MatrixXd I = MatrixXd::Identity(x_.size(),x_.size());
  P_ = (I - K * H_) * P_; 
}
