#pragma once

#include <Eigen/Dense>

#include <iostream>

#include <math.h>

#define DEFAULT_H 0.05

#define LOWER_H_LIMIT_FACTOR 20

#define RMS_ERR_CUTOFF 1.00E-05

#define INITIAL_H_FACTOR 1

class NonlinearSystem {
  double h;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  NonlinearSystem();

  Eigen::VectorXd propagateRK4(double tf, Eigen::VectorXd x0);

  Eigen::VectorXd propagateRK4_adaptive(double tf, Eigen::VectorXd x0);

  void setStepSize(double _h) { h = _h; };

  virtual Eigen::VectorXd f(Eigen::VectorXd x) = 0;
};
