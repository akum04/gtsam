#pragma once

#include <cmath>
#include <iostream>
#include <ostream>

#include <gtsam/geometry/Point3.h>

#include <Eigen/Dense>
#include "math.h"

namespace gtsam {

class principalAxesFrame {
  friend std::ostream& operator<<(std::ostream& out,
                                  const principalAxesFrame& p) {
    p.write(out);
    return out;
  }
  // position	from the target frame	to the	principal frame
  Vector3 _r;
  // quaternion	from	the	target frame to the	principal frame
  Vector4 _q;
  // parameter attitude	error - Modified Rodrigues Parameter
  Vector3 _a;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  static const int dim = 6;
  static const char* name() { return "principalAxesFrame"; }
  Matrix6 _sqrtinf;

  principalAxesFrame() {
    _r << 0.0, 0.0, 0.0;
    _a << 0.0, 0.0, 0.0;
    _q << 0.0, 0.0, 0.0, 1.0;
  }

  principalAxesFrame(Vector3 r) { _r = r; }

  principalAxesFrame(Vector3 r, Vector4 q) {
    _r = r;
    _q = q;
  }

  Vector6 vector() const {
    Vector6 x;
    x << _r, _a;
    return x;
  }

  void set(const Vector6& v) {
    _r = v.block<3, 1>(0, 0);
    _a = v.block<3, 1>(3, 0);
  }

  Vector6 x() {
    Vector6 x;
    x << _r, _a;
    return x;
  }

  Vector4 q() { return _q; }

  Vector4 mrp2quaternion(Vector3 mrp) const {
    Vector4 dq;
    dq << 8 * mrp / (16 + mrp.squaredNorm()),
        (16 - mrp.squaredNorm()) / (16 + mrp.squaredNorm());
    return dq;
  }

  Vector4 addQuaternionError(Vector3 mrp, Vector4 qref) const {
    Vector4 qnew, dq;
    dq = mrp2quaternion(mrp);

    qnew = quaternionMultiplication(dq, qref);
    return qnew;
  }

  principalAxesFrame exmap(const Vector6& Delta) const {
    principalAxesFrame res = *this;
    res._r += Delta.block<3, 1>(0, 0);
    res._a += Delta.block<3, 1>(3, 0);
    return res;
  }

  principalAxesFrame exmap_reset(const Vector6& Delta) {
    principalAxesFrame res = *this;

    res._r += Delta.block<3, 1>(0, 0);
    res._a += Delta.block<3, 1>(3, 0);

    res.write();

    // reset	step
    res._q = addQuaternionError(res._a, res._q);
    res._a = Vector3::Zero();

    printf("inertial	reset\n");

    return res;
  }

  void write(std::ostream& out = std::cout) const {
    out << "	" << _r.transpose();
    out << "	" << _q(0) << "	" << _q(1) << "	" << _q(2) << "	" << _q(3);
    out << "	" << _a.transpose();
    out << std::endl;
  }

  Vector4 quaternionMultiplication(Vector4 q1, Vector4 q2) const {
    // q1	\mult	q2
    Matrix4 qm;
    Vector4 result;
    qm << q1(3), q1(2), -q1(1), q1(0), -q1(2), q1(3), q1(0), q1(1), q1(1),
        -q1(0), q1(3), q1(2), -q1(0), -q1(1), -q1(2), q1(3);

    result = qm * q2;
    result /= result.norm();

    return result;
  }

  Matrix3 rotationMatrix(Vector4 q) const {
    Matrix3 rot;

    rot(0, 0) = q(0) * q(0) - q(1) * q(1) - q(2) * q(2) + q(3) * q(3);
    rot(0, 1) = 2 * (q(0) * q(1) + q(2) * q(3));
    rot(0, 2) = 2 * (q(0) * q(2) - q(1) * q(3));

    rot(1, 0) = 2 * (q(0) * q(1) - q(2) * q(3));
    rot(1, 1) = -q(0) * q(0) + q(1) * q(1) - q(2) * q(2) + q(3) * q(3);
    rot(1, 2) = 2 * (q(2) * q(1) + q(0) * q(3));

    rot(2, 0) = 2 * (q(0) * q(2) + q(1) * q(3));
    rot(2, 1) = 2 * (q(2) * q(1) - q(0) * q(3));
    rot(2, 2) = -q(0) * q(0) - q(1) * q(1) + q(2) * q(2) + q(3) * q(3);

    //	std::cout	<<	"q2rot:	"	<<	q	<<	rot
    //<< std::endl;
    return rot;
  }

  Point3 toPrincipalFrame(const Point3& p_m) const {
    Matrix3 R = rotationMatrix(addQuaternionError(_a, _q));
    Vector3 vecBody = R * (p_m.vector() - _r);
    Point3 p_c(vecBody);

    return p_c;
  }

  Point3 fromPrincipalFrame(const Point3& p_m) const {
    Matrix3 R = rotationMatrix(addQuaternionError(_a, _q));
    Vector3 vecBody = R.transpose() * p_m.vector() + _r;
    Point3 p_c(vecBody);

    return p_c;
  }
};

}  // namespace gtsam