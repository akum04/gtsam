#pragma once

#include <Eigen/Dense>
#include <iostream>
#include "InertiaRatios.h"
#include "NonlinearSystem.h"

using namespace Eigen;

class RigidBodyDynamics : public NonlinearSystem {
  gtsam::InertiaRatios _ir;
  Vector4d _qref;
  Vector3d _r, _v, _a, _w;
  Matrix<double, 6, 6> _Q;
  double _sigma_v, _sigma_w;

  Matrix3d crossProductMat(Vector3d vec);

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  RigidBodyDynamics(gtsam::InertiaRatios ir, double sigma_v, double sigma_w);
  void setMassProperties(gtsam::InertiaRatios ir);
  void setCovProperties(double sigma_v, double sigma_w);
  VectorXd f(VectorXd x);
  void setState(VectorXd x, Vector4d q);
  void setState(VectorXd x);
  void reset_qref();
  Vector4d qref() const { return _qref; };
  Vector4d qTotal() const;
  VectorXd symmMat2Vec(Matrix<double, 12, 12> M);
  Matrix<double, 12, 12> vec2symmMat(VectorXd v);
  Vector4d quaternionFromRot(Matrix3d &R) const;
  Vector4d mrp2quaternion(Vector3d mrp) const;
  Vector3d quaternion2mrp(Vector4d q) const;
  Vector4d addQuaternionError(Vector3d &mrp, Vector4d &qref) const;
  Vector4d quaternionMultiplication(Vector4d &q1, Vector4d &q2) const;
  Vector4d quaternionDivision(Vector4d &q1, Vector4d &q2) const;
  Vector3d diffQuaternion(Vector4d &q, Vector4d &qprev, double dt) const;
  Matrix3d rotationMatrix(Vector4d &q) const;
  Matrix3d getJ() const;
  gtsam::InertiaRatios getIR() const;
  void setIR(gtsam::InertiaRatios ir);
  MatrixXd getBw() const;
  double getSigmaV() const;
  double getSigmaW() const;

  VectorXd x() const;
};
