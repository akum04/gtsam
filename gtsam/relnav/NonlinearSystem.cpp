#include "NonlinearSystem.h"

namespace gtsam {

NonlinearSystem::NonlinearSystem() { h_ = DEFAULT_H; }

// RK4 propagation of the vector
Vector NonlinearSystem::propagateRK4(double tf, Vector x0) {
  Vector k1, k2, k3, k4, x;
  double t = 0;
  double dt;
  bool done = false;
  x = x0;

  while (!done) {
    if (tf - h_ - t > 0) {
      dt = h_;
    } else {
      dt = tf - t;
      done = true;
    }

    k1 = dt * this->f(x);
    k2 = dt * this->f(x + 0.5 * k1);
    k3 = dt * this->f(x + 0.5 * k2);
    k4 = dt * this->f(x + k3);
    x = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    t += dt;
  }

  return x;
}

// RK4 adaptive propagation
Vector NonlinearSystem::propagateRK4_adaptive(double tf, Vector x0) {
  bool done = false;
  double h_starting = this->h_;

  this->h_ = tf / INITIAL_H_FACTOR;

  Vector newX, errX, currX;
  currX = this->propagateRK4(tf, x0);

  while (!done) {
    // new	h_	step	size
    this->h_ = this->h_ / 2;

    // try	the	new	step	size
    newX = this->propagateRK4(tf, x0);

    // compute	rms	error
    errX = newX - currX;
    double rms_err = sqrt(errX.squaredNorm() / errX.size());

    // check	rms_error	or	if	h_	is	too	small
    // that it will take too long
    if (rms_err < RMS_ERR_CUTOFF || this->h_ <= (tf / LOWER_H_LIMIT_FACTOR)) {
      done = true;
      if (this->h_ <= (tf / LOWER_H_LIMIT_FACTOR)) {
      }
    } else {
      currX = newX;
    }
  }

  this->h_ = h_starting;
  return newX;
}

}  // namespace gtsam