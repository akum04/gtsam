#pragma once

#include <iostream>

#include <gtsam/config.h>

#include <Eigen/Dense>

#include <gtsam/geometry/Pose3.h>
#include "InertiaRatios.h"
#include "NonlinearSystem.h"

namespace gtsam {

/**
 * Rigid Body Dynamics with inertia ratios, quaternion,
 * @addtogroup geometry
 * \nosubgrouping
 */

class RigidBody : public NonlinearSystem {
 private:
  InertiaRatios ir_;
  Pose3 pose_;
  Vector3 v_, w_;

 protected:
  Matrix3 SkewSymmetric(Vector3 q) {
    Matrix3 M = Matrix3::Zero();
    M(0, 1) = -q(2);
    M(0, 2) = q(1);
    M(1, 0) = q(2);
    M(1, 2) = -q(0);
    M(2, 0) = -q(1);
    M(2, 1) = q(0);

    return M;
  }

  void RigidBody::setMassProperties(InertiaRatios ir) { ir_ = ir; }

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  RigidBody(InertiaRatios ir);

  void SetState(Vector x, Vector4 q);
  void SetState(Vector x);
  // void reset_qref();
  // Vector4 qTotal() const;

  Vector f(Vector x);
  Vector4 qref() const { return qref_; };
  Vector symmMat2Vec(Eigen::Matrix<double, 12, 12> M);
  Eigen::Matrix<double, 12, 12> vec2symmMat(Vector v);
  // Vector4 quaternionFromRot(Matrix3 &R) const;
  // Vector4 mrp2quaternion(Vector3 mrp) const;
  // Vector3 quaternion2mrp(Vector4 q) const;
  // Vector4 addQuaternionError(Vector3 &mrp, Vector4 &qref) const;
  // Vector4 quaternionMultiplication(Vector4 &q1, Vector4 &q2) const;
  // Vector4 quaternionDivision(Vector4 &q1, Vector4 &q2) const;
  // Vector3 diffQuaternion(Vector4 &q, Vector4 &qprev, double dt) const;
  // Matrix3 rotationMatrix(Vector4 &q) const;
  Matrix3 getJ() const;
  InertiaRatios getIR() const;
  void setIR(InertiaRatios ir);
  Eigen::MatrixXd getBw() const;
  double getSigmaV() const;
  double getSigmaW() const;

  Vector x() const;
  Vector Xq() const;
};

}  // namespace gtsam