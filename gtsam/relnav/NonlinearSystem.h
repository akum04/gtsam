#pragma once

#include <Eigen/Dense>

#include <iostream>

#include <math.h>

#include <gtsam/base/VectorSpace.h>

#define DEFAULT_H 0.05

#define LOWER_H_LIMIT_FACTOR 20

#define RMS_ERR_CUTOFF 1.00E-05

#define INITIAL_H_FACTOR 1

namespace gtsam {

class NonlinearSystem {
  double h;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  NonlinearSystem();

  Vector propagateRK4(double tf, Vector x0);

  Vector propagateRK4_adaptive(double tf, Vector x0);

  void setStepSize(double _h) { h = _h; };

  virtual Vector f(Vector x) = 0;
};

}  // namespace gtsam