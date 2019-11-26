#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>

#include <gtsam/base/concepts.h>
#include <gtsam/geometry/concepts.h>

#include <gtsam/base/VectorSpace.h>

namespace gtsam {

class KinematicPose3 {
 private:
  Vector4 qref_;
  Vector3 r_;
  Vector3 a_;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  static const int dim = 6;
  static const char* name() { return "KinematicPose3"; }

  KinematicPose3() {
    qref_ << 0.0, 0.0, 0.0, 1.0;
    a_ << 0.0, 0.0, 0.0;
    r_ << 0.0, 0.0, 0.0;
  }

  KinematicPose3(const Matrix& hm) {
    // Convert	matrix	to	R,T
    Matrix4 HM = hm / hm(3, 3);  //	enforce	T(3,3)=1
    Matrix3 R = HM.topLeftCorner(3, 3);
    Vector3 r_ = HM.col(3).head(3);

    // compute	quaternion
    qref_ = QuatFromRot(R);
    a_ = Vector3::Zero();
  }

  Vector6 x() const {
    Vector6 x;
    x.segment<3>(0) = r_;
    x.segment<3>(3) = a_;
    return x;
  }

  KinematicPose3(const Vector6& x, const Vector4& q)
      : r_(x.segment<3>(0)), a_(x.segment<3>(3)), qref_(q / q.norm()) {}

  KinematicPose3(const Vector6& x)
      : r_(x.segment<3>(0)), a_(x.segment<3>(3)), qref_(0.0, 0.0, 0.0, 1.0) {}

  Vector4 mrp2quaternion(Vector3 mrp) const {
    Vector4 dq;
    dq << 8 * mrp / (16 + mrp.transpose() * mrp),
        (16 - mrp.transpose() * mrp) / (16 + mrp.transpose() * mrp);
    dq /= dq.norm();
    return dq;
  }

  Vector3 quaternion2mrp(Vector4 q) const {
    Vector3 mrp;
    if (q(3) < 0) {
      q = -q;
    }

    mrp << 4 * q(0) / (1 + q(3)), 4 * q(1) / (1 + q(3)), 4 * q(2) / (1 + q(3));
    return mrp;
  }

  Vector4 addQuaternionError(Vector3& mrp, Vector4& qref) const {
    Vector4 qnew, dq;
    dq = mrp2quaternion(mrp);

    qnew = quaternionMultiplication(dq, qref);

    return qnew;
  }

  Vector4 quaternionMultiplication(Vector4& q1, Vector4& q2) const {
    // q1	\mult	q2
    Matrix4 qm;
    Vector4 result;
    qm << q1(3), q1(2), -q1(1), q1(0), -q1(2), q1(3), q1(0), q1(1), q1(1),
        -q1(0), q1(3), q1(2), -q1(0), -q1(1), -q1(2), q1(3);

    result = qm * q2;
    result /= result.norm();

    return result;
  }

  Vector4 quaternionDivision(Vector4& q1, Vector4& q2) const {
    Vector4 q2inv;

    q2inv << -q2(0), -q2(1), -q2(2), q2(3);

    Vector4 result = quaternionMultiplication(q1, q2inv);
    return result;
  }

  Matrix3 RotationMatrix(Vector4& q) const {
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

    return rot;
  }

  Vector4 QuatFromRot(Matrix3& R) const {
    Vector4 q;
    double div1, div2, div3, div4;

    double numerical_limit = 1.0e-4;

    if (std::abs(R.determinant() - 1) > numerical_limit) {
      std::cerr
          << "R	does	not	have	a	determinant	of	+1"
          << std::endl;
    } else {
      div1 = 0.5 * sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
      div2 = 0.5 * sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
      div3 = 0.5 * sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
      div4 = 0.5 * sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));

      // if	(div1	>	div2	&&	div1	>
      // div3	&&	div1	>	div4)	{
      if (fabs(div1) > numerical_limit) {
        q(3) = div1;
        q(0) = 0.25 * (R(1, 2) - R(2, 1)) / q(3);
        q(1) = 0.25 * (R(2, 0) - R(0, 2)) / q(3);
        q(2) = 0.25 * (R(0, 1) - R(1, 0)) / q(3);
      } else if (fabs(div2) > numerical_limit) {
        //}	else	if	(div2	>	div1	&&	div2
        //>	div3	&&	div2	>	div4)	{
        q(0) = div2;
        q(1) = 0.25 * (R(0, 1) + R(1, 0)) / q(0);
        q(2) = 0.25 * (R(0, 2) + R(2, 0)) / q(0);
        q(3) = 0.25 * (R(1, 2) + R(2, 1)) / q(0);
      } else if (fabs(div3) > numerical_limit) {
        //}	else	if	(div3	>	div1	&&	div3
        //>	div2	&&	div3	>	div4)	{
        q(2) = div3;
        q(0) = 0.25 * (R(0, 2) + R(2, 0)) / q(2);
        q(1) = 0.25 * (R(1, 2) + R(2, 1)) / q(2);
        q(3) = 0.25 * (R(0, 1) - R(1, 0)) / q(2);
        //}	else	{
      } else if (fabs(div4) > numerical_limit) {
        q(1) = div4;
        q(0) = 0.25 * (R(0, 1) + R(1, 0)) / q(1);
        q(2) = 0.25 * (R(1, 2) + R(2, 1)) / q(1);
        q(3) = 0.25 * (R(2, 0) - R(0, 2)) / q(1);
      } else {
        std::cerr << "QuatFromRot	didn't	convert:	[" << div1
                  << ",	" << div2 << ",	" << div3 << ",	" << div4 << std::endl;
        std::cerr << "Rotation	Matrix:	" << R << std::endl;
      }
    }
    q /= q.norm();

    return q;
  }

  Vector3 r() const { return r_; }
  Vector3 a() const { return a_; }
  Vector4 qref() const { return qref_; }

  void ResetQref() {
    Vector3 a_ = a_;
    Vector4 qref_ = qref_;
    qref_ = addQuaternionError(a_, qref_);
    a_ = Vector3::Zero();
  }

  Vector4 qTotal() const {
    Vector3 a_ = a_;
    Vector4 qref_ = qref_;
    return addQuaternionError(a_, qref_);
  };

  KinematicPose3 Expmap(const Vector6& Delta) {
    KinematicPose3 res = *this;
    res.r_ += Delta.head(3);
    res.a_ += Delta.tail(3);
    return res;
  }

  KinematicPose3 ExpmapReset(const Vector6& Delta) {
    KinematicPose3 res = *this;
    res.r_ += Delta.head(3);
    res.a_ += Delta.tail(3);
    res.ResetQref();
    return res;
  }
  // Vector6 vector() const {
  //   Vector6 tmp;
  //   tmp << r_, a_;
  //   return tmp;
  // }

  // void set(const Vector6& v) {
  //   r_ = v.head(3);
  //   a_ = v.tail(3);
  // }

  void write(std::ostream& out) const {
    out << std::endl << "kinPose3	x:	" << x().transpose() << std::endl;
    out << "kinPose3	qref:	" << qref().transpose() << std::endl;
    out << std::endl;
  }

  /**
   *	Convert	Pose3	to	homogeneous	4x4	transformation	matrix.
   *	The	returned	matrix	is	the	object	coordinate
   *frame in	the	world
   *	coordinate	frame.	In	other	words	it	transforms
   *a point	in	the	object frame	to	the	world	frame.
   *
   *	@return	wTo
   */
  Matrix4 wTo() const {
    /*
    Matrix4d	T;
    Vector4	qtot	=	qTotal();
    T.topLeftCorner(3,3)	=	RotationMatrix(qtot).transpose();
    T.col(3).head(3)	=	r_;
    T.row(3)	<<	0.,	0.,	0.,	1.;
    return	T;
    */
    Vector4 qtot = qTotal();
    Matrix3 R = RotationMatrix(qtot);
    Matrix3 oRw = R;
    Vector3 C = -oRw * r_;
    Matrix4 T;
    T.topLeftCorner(3, 3) = oRw;
    T.col(3).head(3) = C;
    T.row(3) << 0., 0., 0., 1.;
    return T;
  }

  /**
   *	Convert	Pose3	to	homogeneous	4x4	transformation matrix.
   *Avoids	inverting	wTo.
   *	The	returned	matrix	is	the	world	coordinate
   *frame in	the	object
   *	coordinate	frame.	In	other	words	it	transforms
   *a point	in	the	world frame	to	the	object	frame.
   *
   *	@return	oTw
   */
  Matrix4 oTw() const {
    Matrix4 T;
    Vector4 qtot = qTotal();
    T.topLeftCorner(3, 3) = RotationMatrix(qtot).transpose();
    T.col(3).head(3) = r_;
    T.row(3) << 0., 0., 0., 1.;
    return T;
  }

  bool equals(const KinematicPose3& pose, double tol) const {
    return traits<Vector3>::Equals(a_, pose.a_, tol) &&
           traits<Vector3>::Equals(r_, pose.r_, tol);
  }

  void Print(const std::string& str = " ") {
    if (str.size() == 0)
      std::cout << "Kinematic Pose: ";
    else
      std::cout << str << " ";
    std::cout << "kinPose3	x:	" << x().transpose() << std::endl;
    std::cout << "kinPose3	qref:	" << qref().transpose() << std::endl;
  }
};

std::ostream& operator<<(std::ostream& out, const KinematicPose3& p) {
  out << std::endl << "kinPose3	x:	" << p.x().transpose() << std::endl;
  out << "kinPose3	qref:	" << p.qref().transpose() << std::endl;
  out << std::endl;
  return out;
}

}  // namespace gtsam
