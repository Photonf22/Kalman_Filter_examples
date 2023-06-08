#include <iostream>
#include "ukf.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

UKF::UKF() {
  Init();
}

UKF::~UKF() {

}

void UKF::Init() {

}


/**
 * Programming assignment functions: 
 */

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

  // set state dimension
  int n_x = 5;

  // set augmented dimension
  int n_aug = 7;

  // // create example sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_aug <<
    5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
      1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
    2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
    0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
    0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
         0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
         0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;
  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  double delta_t = 0.1; // time diff in sec

  /**
   * Student part begin
   */
  // predict sigma points
    VectorXd state_prt0 =  VectorXd(5);
    VectorXd state_prt1 = VectorXd(5);
    VectorXd state_prt2 = VectorXd(5);
  for(int i = 0 ;i < (2 * n_aug + 1);i++)
  {
    state_prt0 <<  Xsig_aug(0,i),Xsig_aug(1,i),Xsig_aug(2,i),Xsig_aug(3,i),Xsig_aug(4,i);
    
    if(Xsig_aug(4,i) == 0)
    {
      state_prt1<< Xsig_aug(2,i)*cos(Xsig_aug(3,i))*delta_t,Xsig_aug(2,i)*sin(Xsig_aug(3,i))*delta_t,(0),(Xsig_aug(4,i)*delta_t),(0);
      state_prt2<< (0.5)*((delta_t*delta_t)*cos(Xsig_aug(3,i))*Xsig_aug(5,i)),(0.5)*((delta_t*delta_t)*sin(Xsig_aug(3,i))*Xsig_aug(5,i)),(delta_t*Xsig_aug(5,i)),(1/2)*((delta_t*delta_t)*Xsig_aug(6,i)),(delta_t*Xsig_aug(6,i));
    }
    else
    {
      state_prt1<< (Xsig_aug(2,i)/Xsig_aug(4,i)) * (sin(Xsig_aug(3,i) + Xsig_aug(4,i)*delta_t) - sin(Xsig_aug(3,i))),(Xsig_aug(2,i)/Xsig_aug(4,i)) * (-cos(Xsig_aug(3,i) + Xsig_aug(4,i)*delta_t) + cos(Xsig_aug(3,i))),(0),(Xsig_aug(4,i)*delta_t),(0);
      state_prt2<< (0.5)*((delta_t*delta_t)*cos(Xsig_aug(3,i))*Xsig_aug(5,i)),(0.5)*((delta_t*delta_t)*sin(Xsig_aug(3,i))*Xsig_aug(5,i)),(delta_t*Xsig_aug(5,i)),(0.5)*((delta_t*delta_t)*Xsig_aug(6,i)),(delta_t*Xsig_aug(6,i));
    }
 
    Xsig_pred.col(i) =  state_prt0 + state_prt1 + state_prt2;
     
  }
  // print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  // write result
  *Xsig_out = Xsig_pred;
}