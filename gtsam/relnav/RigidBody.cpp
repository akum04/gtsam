/* ----------------------------------------------------------------------------

* Author : Arunkumar Rathinam

 * -------------------------------------------------------------------------- */

/**
 * @file  RigidBody.cpp
 * @brief RigidBody Dynamcis
 */

// #include <gtsam/geometry/Pose3.h>
#include <gtsam/relnav/RigidBody.h>

namespace gtsam {

// RigidBody::RigidBody() {
//   pose_ = Pose3();
//   v_ = Vector3::Zero();
//   w_ = Vector3::Zero();
// }

RigidBody::RigidBody(Pose3 pose, Vector3 v, Vector3 w) {
  pose_ = pose;
  v_ = v;
  w_ = w;
}

/* ************************************************************************ */

void RigidBody::SetPose(Pose3 pose) {
  Rot3 r = pose.rotation();
  Vector3 t = pose.translation();
  pose_ = Pose3(r, t);
}

Pose3 RigidBody::GetPose() { return pose_; }

/* ************************************************************************ */

void RigidBody::SetState(Vector x) {
  if (x.size() == 12) {
    pose_ = Pose3(Rot3::RzRyRx(x.segment<3>(3)), x.segment<3>(0));
    v_ = x.segment<3>(6);
    w_ = x.segment<3>(9);
  }
}

/* ********************** */

// void RigidBody::SetState(Vector x, Vector4 q) {
//   if (x.size() == 12) {
//     pose_ = Pose3(Rot3::Quaternion(q(0), q(1), q(2), q(3)), x.segment<3>(0));
//     v_ = x.segment<3>(6);
//     w_ = x.segment<3>(9);
//   }
// }

/* ********************** */

// void RigidBody::SetState(Vector x, Quaternion q) {
//   if (x.size() == 12) {
//     pose_ = Pose3(q, x.segment<3>(0));
//     v_ = x.segment<3>(6);
//     w_ = x.segment<3>(9);
//   }
// }

/* ************************************************************************ */

Vector RigidBody::f(Vector x) {
  Vector3 dr, dv, da, dw;
  Eigen::Matrix<double, 12, 12> lambda, dLambda;
  Vector vec_dLambda;
  Vector dx(90);

  Vector3 r = x.segment<3>(0);
  Vector3 a = x.segment<3>(3);
  Vector3 v = x.segment<3>(6);
  Vector3 w = x.segment<3>(9);

  Matrix Bw = getBw();
  Matrix3 J = Matrix3::Identity();
  // Matrix3 J = ir_.getJ();

  // Nonlinear	State	Model	\dot	x	=	f(x)

  /*
   *	\mathbf{\dot	r}	=	\mathbf{v}
   */
  dr = v;

  /*
   *	\mathbf{\dot	v}	=	0
   */
  dv = Vector3::Zero();

  /*
   *	\frac{d	\mathbf{a}_p}{dt}	=
   *	\frac{1}{2}\left(\mathbf{[\omega	\times]}	+
   *	\mathbf{\omega}	\cdot	\mathbf{\bar	q}	\right) \mathbf{a}_p
   *+ \frac{2	q_4}{1+q_4}	\mathbf{\omega}
   */
  double c1, c2, c3;
  c1 = 0.5;
  c2 = 0.125 * w.dot(a);
  c3 = 1 - a.dot(a) / 16;
  da = -c1 * w.cross(a) + c2 * a + c3 * w;

  /*
   *	\dot	\mathbf{w}	=	-\mathbf{J}ˆ{-1}	\mathbf{\omega}
   *\times \mathbf{J}	\mathbf{\omega}
   */
  dw = -J.inverse() * w.cross(J * w);
  /*
  // Covariance	Propagation	according	to	Lyapunov
  // function see	Brown	&	Hwang	pg	204

  // Compute	Linear	transition	matrix
  Eigen::Matrix<double, 12, 12> A = Eigen::Matrix<double, 12, 12>::Zero();

  // position	derivative
  A.block<3, 3>(0, 3) = Matrix3::Identity();

  // mrp	kinematics
  A.block<3, 3>(6, 6) =
      -0.5 * SkewSymmetric(w) + w.dot(a) / 8 * Matrix3::Identity();
  A.block<3, 3>(6, 9) = (1 - a.dot(a / 16)) * Matrix3::Identity();

  // angular	velocity	dynamics
  A.block<3, 3>(9, 9) = -J.inverse() * SkewSymmetric(w) * J;

  lambda = vec2symmMat(x.segment<78>(12));
  dLambda = A * lambda + lambda * A.transpose() + Bw * Q_ * Bw.transpose();
  vec_dLambda = symmMat2Vec(dLambda);
*/
  // write	to	dx
  dx.segment<3>(0) = dr;
  dx.segment<3>(3) = da;
  dx.segment<3>(6) = dv;
  dx.segment<3>(9) = dw;
  // dx.segment<78>(12) = vec_dLambda;

  return dx;
}

/* ************************************************************************ */

Vector RigidBody::symmMat2Vec(Eigen::Matrix<double, 12, 12> M) {
  Vector v(78);
  int count = 0;
  for (int row = 0; row < 12; row++) {
    for (int col = row; col < 12; col++) {
      v(count) = M(row, col);
      count++;
    }
  }
  return v;
}

/* ************************************************************************ */

Eigen::Matrix<double, 12, 12> RigidBody::vec2symmMat(Vector v) {
  Eigen::Matrix<double, 12, 12> M = Eigen::Matrix<double, 12, 12>::Zero();
  int count = 0;
  for (int row = 0; row < 12; row++) {
    for (int col = row; col < 12; col++) {
      M(row, col) = v(count);
      M(col, row) = v(count);
      count++;
    }
  }
  return M;
}

/* ************************************************************************ */

// Vector RigidBody::Xq() const {
//   Vector x(13);
//   x.segment<3>(0) = pose_.translation();
//   Quaternion q = pose_.rotation().toQuaternion();
//   x(3) = q.w();
//   x.segment<3>(4) = q.vec();
//   x.segment<3>(7) = v_;
//   x.segment<3>(10) = w_;
//   return x;
// }
/* ************************************************************************ */

Vector RigidBody::x() const {
  Vector x(12);
  x.segment<3>(0) = pose_.translation();
  x.segment<3>(3) = pose_.rotation().xyz();
  x.segment<3>(6) = v_;
  x.segment<3>(9) = w_;
  return x;
}

/* ************************************************************************ */

Matrix RigidBody::getBw() const {
  Eigen::Matrix<double, 12, 6> Bw;
  Bw = Eigen::Matrix<double, 12, 6>::Zero();
  Bw.block<3, 3>(3, 0) = Matrix3::Identity();
  Bw.block<3, 3>(9, 3) = Matrix3::Identity();

  return Bw;
}

/* ************************************************************************ */

// Vector4 RigidBody::qTotal() const {
//   Vector3 a = a_;
//   Vector4 qref;
//   qref << 1.0, 0.0, 0.0, 0.0;
//   // Quaternion qr = Quaternion(qref(0), qref(1), qref(2), qref(3));
//   // Quaternion qa = Rot3::RzRyRx(a(0), a(1), a(2)).toQuaternion();
//   // Quaternion qn = qa * qr;
//   // if (qn.w() * qr.w() + qn.vec().dot(qr.vec()) >= 0) {
//   //   Vector4 Q = (qn.w(), qn.vec());
//   //   return Q;
//   // } else {
//   //   Vector3 qnv = -1.0 * qn.vec();
//   //   double qnv0 = qnv(0);
//   //   double qnw = -qn.w();
//   //   Quaternion qn2 = Quaternion(qnw, qnv0, qnv0, qnv0);
//   //   Vector4 Q = (qn2.w(), qn2.vec());
//   //   return Q;
//   // }
//   return qref;  // addQuaternionError(a_, qref_);
// }

/* ************************************************************************ */

void RigidBody::print(const std::string& s) {
  std::cout << " Rigid body: \n";
  pose_.print("Pose: \n");
  std::cout << "Linear and Angular Velocity: \n"
            << x().segment<3>(6).transpose() << std::endl;
  // ir_.print("Inertia Ratios: \n");
}

/* ************************************************************************* */
bool RigidBody::equals(const RigidBody& rb, double tol = 1e-9) const {
  return pose_.rotation().equals(rb.pose_.rotation(), tol) &&
         traits<Point3>::Equals(pose_.translation(), rb.pose_.translation(),
                                tol) &&
         traits<Point3>::Equals(v_, rb.v_, tol) &&
         traits<Point3>::Equals(w_, rb.w_, tol);
}

RigidBody RigidBody::between(const RigidBody& prev,
                             OptionalJacobian<12, 12> H1 = boost::none,
                             OptionalJacobian<12, 12> H2 = boost::none) const {
  Pose3 pose;
  // Pose3 pose = traits<Pose3>::Between(pose_, prev.pose_, H1<6, 6>(0, 0),
  //                                     H2<6, 6>(0, 0));  // h(x)

  // Pose3 pose =
  //     pose_.Between(prev.pose_, H1<6, 6>(0, 0), H2<6, 6>(0, 0));  // h(x)
  Vector6 v, vw, vwp;
  vw << v_, w_;
  vwp << prev.v_, prev.w_;
  // v = traits<Vector6>::Between(vw, vwp, H1<6, 6>(6, 6), H2<6, 6>(6, 6));
  RigidBody rb(pose, v.segment<3>(0), v.segment<3>(3));
  return rb;
}

RigidBody RigidBody::localCoordinates(const RigidBody& prev) const {
  Pose3 pose;
  // Pose3 pose = traits<Pose3>::Local(pose_, prev.pose_);
  Vector6 v, vw, vwp;
  // vw << v_, w_;
  // vwp << prev.v_, prev.w_;
  // v = traits<Vector6>::Local(vw, vwp);
  RigidBody rb(pose, v.segment<3>(0), v.segment<3>(3));
  return rb;
}

}  // namespace gtsam
