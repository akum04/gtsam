#include "NonlinearSystem.h"

NonlinearSystem::NonlinearSystem() { h = DEFAULT_H; }

// RK4 propagation of the vector
Eigen::VectorXd NonlinearSystem::propagateRK4(double tf, Eigen::VectorXd x0) {
  Eigen::VectorXd k1, k2, k3, k4, x;
  double t = 0;
  double dt;
  bool done = false;
  x = x0;

  while (!done) {
    if (tf - h - t > 0) {
      dt = h;
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
Eigen::VectorXd NonlinearSystem::propagateRK4_adaptive(double tf,
                                                       Eigen::VectorXd x0) {
  bool done = false;
  double h_starting = this->h;

  this->h = tf / INITIAL_H_FACTOR;

  Eigen::VectorXd newX, errX, currX;
  currX = this->propagateRK4(tf, x0);

  while (!done) {
    // new	h	step	size
    this->h = this->h / 2;

    // try	the	new	step	size
    newX = this->propagateRK4(tf, x0);

    // compute	rms	error
    errX = newX - currX;
    double rms_err = sqrt(errX.squaredNorm() / errX.size());

    // check	rms_error	or	if	h	is	too	small
    // that it will take too long
    if (rms_err < RMS_ERR_CUTOFF || this->h <= (tf / LOWER_H_LIMIT_FACTOR)) {
      done = true;
      if (this->h <= (tf / LOWER_H_LIMIT_FACTOR)) {
      }
    } else {
      currX = newX;
    }
  }

  this->h = h_starting;
  return newX;
}
