#pragma once

#include <iostream>

#include <Eigen/Dense>

#include <gtsam/base/concepts.h>
#include <gtsam/config.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/concepts.h>

#include "InertiaRatios.h"
#include "NonlinearSystem.h"

namespace gtsam {

/**
 * Rigid Body Dynamics with inertia ratios, quaternion,
 * @addtogroup geometry
 * \nosubgrouping
 */

/**
 * A copy of Rigidbody file is in archive folder
 *
 */

class RigidBody : public NonlinearSystem {
 private:
  Pose3 pose_;
  Vector3 v_, w_;

 protected:
  Matrix3 SkewSymmetric(Vector3 q) {
    Matrix3 M;
    M << 0.0, -q(2), q(1), q(2), 0.0, -q(0), -q(1), q(0), 0.0;
    return M;
  }

  // void SetMassProperties(InertiaRatios ir) { ir_ = ir; }

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  RigidBody(Pose3 pose, Vector3 v, Vector3 w);

  void SetPose(Pose3 pose);
  Pose3 GetPose();

  void SetState(Vector x);
  // void SetState(Vector x, Vector4 q);
  // void SetState(Vector x, Quaternion q);

  Vector f(Vector x);

  Vector symmMat2Vec(Eigen::Matrix<double, 12, 12> M);
  Eigen::Matrix<double, 12, 12> vec2symmMat(Vector v);
  Eigen::MatrixXd getBw() const;

  Vector x() const;

  void print(const std::string& s);
  bool equals(const RigidBody& rb, double tol) const;

  RigidBody between(const RigidBody& prev, OptionalJacobian<12, 12> H1,
                    OptionalJacobian<12, 12> H2) const;

  RigidBody localCoordinates(const RigidBody& prev) const;
};

}  // namespace gtsam