#pragma once

#include <iostream>

#include <Eigen/Dense>

#include "InertiaRatios.h"
#include "NonlinearSystem.h"

namespace gtsam {

class RigidBodyDynamics : public NonlinearSystem {
 private:
  InertiaRatios ir_;
  Vector4 qref_;
  Vector3 r_, v_, a_, w_;
  Matrix6 Q_;
  double sigma_v_, sigma_w_;

  Matrix3 crossProductMat(Vector3 vec);

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  RigidBodyDynamics(InertiaRatios ir, double sigma_v, double sigma_w);
  void setMassProperties(InertiaRatios ir);
  void setCovProperties(double sigma_v, double sigma_w);
  Vector f(Vector x);
  void setState(Vector x, Vector4 q);
  void setState(Vector x);
  void reset_qref();
  Vector4 qref() const { return qref_; };
  Vector4 qTotal() const;
  Vector symmMat2Vec(Eigen::Matrix<double, 12, 12> M);
  Eigen::Matrix<double, 12, 12> vec2symmMat(Vector v);
  Vector4 quaternionFromRot(Matrix3 &R) const;
  Vector4 mrp2quaternion(Vector3 mrp) const;
  Vector3 quaternion2mrp(Vector4 q) const;
  Vector4 addQuaternionError(Vector3 &mrp, Vector4 &qref) const;
  Vector4 quaternionMultiplication(Vector4 &q1, Vector4 &q2) const;
  Vector4 quaternionDivision(Vector4 &q1, Vector4 &q2) const;
  Vector3 diffQuaternion(Vector4 &q, Vector4 &qprev, double dt) const;
  Matrix3 rotationMatrix(Vector4 &q) const;
  Matrix3 getJ() const;
  InertiaRatios getIR() const;
  void setIR(InertiaRatios ir);
  Eigen::MatrixXd getBw() const;
  double getSigmaV() const;
  double getSigmaW() const;

  Vector x() const;
};

}  // namespace gtsam