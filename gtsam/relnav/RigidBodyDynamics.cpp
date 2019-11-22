#include "RigidBodyDynamics.h"

namespace gtsam {

RigidBodyDynamics::RigidBodyDynamics(InertiaRatios ir, double sigma_v,
                                     double sigma_w) {
  ir_ = ir;
  qref_ << 0, 0, 0, 1;
  r_ = Vector3::Zero();
  v_ = Vector3::Zero();
  a_ = Vector3::Zero();
  w_ = Vector3::Zero();

  setMassProperties(ir);
  setCovProperties(sigma_v, sigma_w);
}

void RigidBodyDynamics::setMassProperties(InertiaRatios ir) { ir_ = ir; }

void RigidBodyDynamics::setCovProperties(double sigma_v, double sigma_w) {
  sigma_v_ = sigma_v;
  sigma_w_ = sigma_w;
  Q_ = Matrix6::Zero();
  Q_.block<3, 3>(0, 0) = sigma_v_ * sigma_v_ * Matrix3::Identity();
  Q_.block<3, 3>(3, 3) = sigma_w_ * sigma_w_ * Matrix3::Identity();
}

void RigidBodyDynamics::reset_qref() {
  Vector3 a_ = a_;
  Vector4 qref_ = qref_;
  qref_ = addQuaternionError(a_, qref_);
  a_ = Vector3::Zero();
}

Vector4 RigidBodyDynamics::qTotal() const {
  Vector3 a_ = a_;
  Vector4 qref_ = qref_;
  return addQuaternionError(a_, qref_);
};

Vector RigidBodyDynamics::f(Vector x) {
  Vector3 dr, dv, da, dw;
  Eigen::Matrix<double, 12, 12> lambda, dLambda;
  Vector vec_dLambda;
  Vector dx(90);

  Vector3 r = x.segment<3>(0);
  Vector3 v = x.segment<3>(3);
  Vector3 a = x.segment<3>(6);
  Vector3 w = x.segment<3>(9);

  Matrix Bw = getBw();
  Matrix3 J = ir_.getJ();

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

  // Covariance	Propagation	according	to	Lyapunov	function
  // see	Brown	&	Hwang	pg	204

  // Compute	Linear	transition	matrix
  Eigen::Matrix<double, 12, 12> A = Eigen::Matrix<double, 12, 12>::Zero();

  // position	derivative
  A.block<3, 3>(0, 3) = Matrix3::Identity();

  // mrp	kinematics
  A.block<3, 3>(6, 6) =
      -0.5 * crossProductMat(w) + w.dot(a) / 8 * Matrix3::Identity();
  A.block<3, 3>(6, 9) = (1 - a.dot(a / 16)) * Matrix3::Identity();

  // angular	velocity	dynamics
  A.block<3, 3>(9, 9) = -J.inverse() * crossProductMat(w) * J;

  lambda = vec2symmMat(x.segment<78>(12));
  dLambda = A * lambda + lambda * A.transpose() + Bw * Q_ * Bw.transpose();
  vec_dLambda = symmMat2Vec(dLambda);
  // write	to	dx
  dx.segment<3>(0) = dr;
  dx.segment<3>(3) = dv;
  dx.segment<3>(6) = da;
  dx.segment<3>(9) = dw;
  dx.segment<78>(12) = vec_dLambda;

  return dx;
}

Matrix3 RigidBodyDynamics::crossProductMat(Vector3 vec) {
  Matrix3 M = Matrix3::Zero();
  M(0, 1) = -vec(2);
  M(0, 2) = vec(1);
  M(1, 0) = vec(2);
  M(1, 2) = -vec(0);
  M(2, 0) = -vec(1);
  M(2, 1) = vec(0);

  return M;
}

Vector RigidBodyDynamics::symmMat2Vec(Eigen::Matrix<double, 12, 12> M) {
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

Eigen::Matrix<double, 12, 12> RigidBodyDynamics::vec2symmMat(Vector v) {
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

Vector RigidBodyDynamics::x() const {
  Vector x(12);
  x.segment<3>(0) = r_;
  x.segment<3>(3) = v_;
  x.segment<3>(6) = a_;
  x.segment<3>(9) = w_;
  return x;
}

void RigidBodyDynamics::setState(Vector x, Vector4 q) {
  r_ = x.segment<3>(0);
  v_ = x.segment<3>(3);
  a_ = x.segment<3>(6);
  w_ = x.segment<3>(9);
  qref_ = q / q.norm();
}

void RigidBodyDynamics::setState(Vector x) {
  r_ = x.segment<3>(0);
  v_ = x.segment<3>(3);
  a_ = x.segment<3>(6);
  w_ = x.segment<3>(9);
}

Vector4 RigidBodyDynamics::mrp2quaternion(Vector3 mrp) const {
  Vector4 dq;
  dq << 8 * mrp / (16 + mrp.transpose() * mrp),
      (16 - mrp.transpose() * mrp) / (16 + mrp.transpose() * mrp);
  dq /= dq.norm();

  return dq;
}

Vector3 RigidBodyDynamics::quaternion2mrp(Vector4 q) const {
  Vector3 mrp;
  if (q(3) < 0) {
    q = -q;
  }

  mrp << 4 * q(0) / (1 + q(3)), 4 * q(1) / (1 + q(3)), 4 * q(2) / (1 + q(3));
  return mrp;
}

Vector4 RigidBodyDynamics::addQuaternionError(Vector3 &mrp,
                                              Vector4 &qref) const {
  Vector4 qnew, dq;
  dq = mrp2quaternion(mrp);

  Vector4 qnew1 = quaternionMultiplication(dq, qref);

  if (qnew1.dot(qref) >= 0) {
    return qnew1;
  } else {
    Vector4 qnew2 = -1 * qnew1;
    return qnew2;
  }
}

Vector4 RigidBodyDynamics::quaternionMultiplication(Vector4 &q1,
                                                    Vector4 &q2) const {
  // q1	\mult	q2
  Matrix4 qm;
  Vector4 result;
  qm << q1(3), q1(2), -q1(1), q1(0), -q1(2), q1(3), q1(0), q1(1), q1(1), -q1(0),
      q1(3), q1(2), -q1(0), -q1(1), -q1(2), q1(3);

  result = qm * q2;
  result /= result.norm();

  return result;
}

Vector4 RigidBodyDynamics::quaternionDivision(Vector4 &q1, Vector4 &q2) const {
  Vector4 q2inv;

  q2inv << -q2(0), -q2(1), -q2(2), q2(3);

  Vector4 result = quaternionMultiplication(q1, q2inv);
  return result;
}

Vector3 RigidBodyDynamics::diffQuaternion(Vector4 &q, Vector4 &qprev,
                                          double dt) const {
  Vector4 dq = (q - qprev) / dt;
  Matrix4 M;

  M << qprev(3), qprev(2), -qprev(1), -qprev(0), -qprev(2), qprev(3), qprev(0),
      -qprev(1), qprev(1), -qprev(0), qprev(3), -qprev(2), qprev(0), qprev(1),
      qprev(2), qprev(3);

  Vector4 wp = 2 * M * dq;
  Vector3 w = wp.head(3);

  return w;
}

Matrix3 RigidBodyDynamics::rotationMatrix(Vector4 &q) const {
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

Vector4 RigidBodyDynamics::quaternionFromRot(Matrix3 &R) const {
  Vector4 q;
  double div1, div2, div3, div4;

  double numerical_limit = 1.0e-4;

  if (std::abs(R.determinant() - 1) > numerical_limit) {
    std::cerr << "R does not have a determinant of + 1" << std::endl;
  } else {
    div1 = 0.5 * sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
    div2 = 0.5 * sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
    div3 = 0.5 * sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
    div4 = 0.5 * sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));

    // if	(div1	>	div2	&&	div1	>	div3	&&
    // div1 > div4)
    // {
    if (fabs(div1) > numerical_limit) {
      q(3) = div1;
      q(0) = 0.25 * (R(1, 2) - R(2, 1)) / q(3);
      q(1) = 0.25 * (R(2, 0) - R(0, 2)) / q(3);
      q(2) = 0.25 * (R(0, 1) - R(1, 0)) / q(3);
    } else if (fabs(div2) > numerical_limit) {
      //}	else	if	(div2	>	div1	&&
      // div2	>	div3	&&	div2	>	div4)	{
      q(0) = div2;
      q(1) = 0.25 * (R(0, 1) + R(1, 0)) / q(0);
      q(2) = 0.25 * (R(0, 2) + R(2, 0)) / q(0);
      q(3) = 0.25 * (R(1, 2) + R(2, 1)) / q(0);
    } else if (fabs(div3) > numerical_limit) {
      //}	else	if	(div3	>	div1	&&
      // div3	>	div2	&&	div3	>	div4)	{
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
      std::cerr << "quaternionFromRot didn't convert: [" << div1 << "," << div2
                << "," << div3 << "," << div4 << std::endl;
      std::cerr << "Rotation Matrix :" << R << std::endl;
    }
  }
  q /= q.norm();

  return q;
}

Matrix RigidBodyDynamics::getBw() const {
  Eigen::Matrix<double, 12, 6> Bw;
  Bw = Eigen::Matrix<double, 12, 6>::Zero();
  Bw.block<3, 3>(3, 0) = Matrix3::Identity();
  Bw.block<3, 3>(9, 3) = Matrix3::Identity();

  return Bw;
}

Matrix3 RigidBodyDynamics::getJ() const { return ir_.getJ(); }

InertiaRatios RigidBodyDynamics::getIR() const { return ir_; }

void RigidBodyDynamics::setIR(InertiaRatios ir) { ir_ = ir; }

double RigidBodyDynamics::getSigmaV() const { return sigma_v_; }

double RigidBodyDynamics::getSigmaW() const { return sigma_w_; }

}  // namespace gtsam