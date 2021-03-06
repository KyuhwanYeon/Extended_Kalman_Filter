#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict()
{
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const Eigen::VectorXd &z, Eigen::MatrixXd &Hj_, Eigen::MatrixXd &R_)
{
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  double rho = sqrt(px * px + py * py);
  rho = (rho < 0.01) ? 0.01 : rho; // to avoid divide by 0
  double phi = atan2(py, px);
  double rho_dot = (px * vx + py * vy) / rho;
  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot; // z_pred = h(x)

  VectorXd y = z - z_pred;

  // In this project, z(1) is sometimes over the pi. Actually in the real world, rho cannot be over the pi.
  // Below saturation code is just dedicated for this project. 
  while (y(1) > M_PI)
  {
    printf("z: %lf, z_pred: %lf\n", z(1), z_pred(1));
    y(1) -= 2. * M_PI;
  }
  while (y(1) < -M_PI)
  {
    printf("z: %lf, z_pred: %lf\n", z(1), z_pred(1));
    y(1) += 2. * M_PI;
  }
  
  MatrixXd Ht = Hj_.transpose();
  MatrixXd S = Hj_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj_) * P_;
}
