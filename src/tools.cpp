#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  
  VectorXd rmse(4);
  rmse = VectorXd::Zero(4);
  
  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
    cout << "Error: Size doesn't match" << endl;
    return rmse;
  }
 
  for(int i=0;i<estimations.size();i++){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array()*residual.array();
    rmse += residual;
  }
  
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt(); 
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  
  MatrixXd Hj(3,4);
  Hj = MatrixXd::Zero(3,4); 
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float d2 = px*px + py*py;
  float d = sqrt(d2);
  float d3 = d2*d;
  
  if(fabs(d)<0.0001){
    cout<<"Error: Division by Zero"<<endl;
    return Hj;
  }
  Hj << (px/d),(py/d), 0,0,
      (-py/d2),(px/d2),0,0,
     py*(vx*py-vy*px)/d3,px*(vy*px-vx*py)/d3,px/d,py/d;
  
  return Hj;
}
