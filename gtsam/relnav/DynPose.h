#pragma once

#include <iostream>

#include <Eigen/Dense>

#include "RigidBody.h"

// using namespace Eigen;
namespace gtsam {

class DynPose {
 private:
  RigidBody rbd_;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  //	assignment	operator	and	copy	constructor
  // implicitly created,	which	is	ok
  static const int dim = 12;
  static const char* name() { return "DynamicPose"; }

  /* ************************************************************************ */

  DynPose(InertiaRatios ir) : rbd_(ir) {}

  // copy	constructor
  DynPose(const DynPose& dp) : rbd_(dp.rbd_.getIR()) {
    rbd_.SetState(dp.rbd_.Xq());
  }

  /* ************************************************************************ */

  DynPose& operator=(const DynPose& dp) {
    rbd_ = RigidBody(dp.rbd_.getIR());
    rbd_.SetPose(dp.rbd_.GetPose());
    return *this;
  }

  /* ************************************************************************ */

  DynPose(Vector x, InertiaRatios ir) : rbd_(ir) { rbd_.SetState(x); }

  DynPose(Vector x, Vector4 qref, InertiaRatios ir) : rbd_(ir) {
    rbd_.SetState(x, qref);
  }

  /* ************************************************************************ */

  DynPose(const Matrix4& hm, bool initVelocities, double dt, InertiaRatios ir)
      : rbd_(ir) {
    Vector12 x;
    Vector3 r;
    Vector3 v;
    Vector3 a;
    // Vector4 q;
    Vector3 w;

    // Convert	matrix	to	R,T
    Matrix4 HM = hm / hm(3, 3);  //	enforce	T(3,3)=1
    Matrix3 R = HM.topLeftCorner(3, 3);
    Vector3 T = HM.col(3).head(3);

    // compute	quaternion
    Rot3 q = Rot3(R(0, 0), R(0, 1), R(0, 2), R(1, 0), R(1, 1), R(1, 2), R(2, 0),
                  R(2, 1), R(2, 2));
    Quaternion qn = q.toQuaternion();
    // q = rbd_.quaternionFromRot(R);
    a = Vector3::Zero();
    Pose3 p = Pose3(q, a);

    if (initVelocities) {
      // differentiate	linear	velocity
      v = T / dt;

      /*	Differentiate	quaternion:
       *	dq/dt	=	(q[k]	-	q[k-1])/dt	= 0.5
       *O(w[k-1])	q[k-1]
       *	where	O(w[k-1])	is	an	orthonormal
       *quaternion	mult	matrix	for	[w1;	w2;
       *w3;	0]	(i.e.	quaternionMultiplication()) set	q[k-1]	=
       *[0;0;0;1]	by	definition	(from	a refererence frame)
       *and	solve	for	w[k-1]	gives
       *	w[k-1]	=	2	[q1[k];	q2[k];	q3[k]]	/	dt
       */
      w = 2 * qn.vec() / dt;
    } else {
      v = Vector3::Zero();
      w = Vector3::Zero();
    }

    x.block<3, 1>(0, 0) = v;
    x.block<3, 1>(3, 0) = a;
    x.block<3, 1>(6, 0) = w;
    // x.block<3, 1>(9, 0) = w;
    rbd_.setState(x, q);
  }

  Vector x() { return rbd_.x(); };

  Vector4 q() {
    Pose3 p = rbd_.GetPose();
    q = p.rotation().toQuaternion();
    x = q.w();
    x.segment<3>(1) = q.vec();
    return x;
  };

  /* ************************************************************************ */

  Vector4 qTotal() const { return rbd_.qTotal(); };

  /* ************************************************************************ */

  DynPose Expmap(const Vector12& Delta) const {
    DynPose res = *this;
    res.rbd_.setState(res.rbd_.x() + Delta);
    return res;
  }

  /* ************************************************************************ */

  DynPose Expmap_reset(const Vector12& Delta) {
    DynPose res = *this;
    res.rbd_.setState(res.rbd_.x() + Delta);
    res.rbd_.reset_qref();

    /*	We	should	REALLY,	REALLY	update	Factor::_sqrtinf	at
     *this location in	the	code.	However	it	is	a	const
     *variable,	and	there	is	no
     *way to	do	callbacks	to	the	Factor	class.
     *So	I	will	leave	this	for	future
     *	work.	Now,	the	value	that	is	created	on
     *initialization is the final version,	even	after	many
     *relinearizations.
     */

    // NOTE	-	THIS	IS	NOW	DONE	in
    // NodeExpmapT->apply_reset();

    return res;
  }

  void set(const Vector& v) { rbd_.setState(v); }

  void set_qref(const Vector4& qset) { rbd_.setState(rbd_.x(), qset); }

  void rezero() {
    Vector x = Vector::Zero(12);
    Vector4 q;
    q << 0, 0, 0, 1;
    rbd_.setState(x, q);
  }

  /* ************************************************************************ */

  DynPose propagate(double dt, InertiaRatios& ir) {
    Vector x0 = Vector::Zero(90);
    x0.head(12) = rbd_.x();
    rbd_.setIR(ir);
    //	std::cout	<<	"x0:	"	<<	x0.head(12).transpose()
    //<< std::endl;
    Vector newX = rbd_.propagateRK4_adaptive(dt, x0).head(12);

    //	std::cout	<<	"dt:	"	<<	dt	<<
    // std::endl; 	std::cout	<<	"newX:	"	<<
    // newX.transpose()
    //<<	std::endl;

    DynPose xNew(newX, this->rbd_.qref(), this->rbd_.getIR(),
                 this->rbd_.getSigmaV(), this->rbd_.getSigmaW());
    xNew.Expmap(Vector12::Zero());
    return xNew;
  }

  /* ************************************************************************ */

  // Vector3 getMRPdifference(Vector4 qtot1, Vector4 qtot2) {
  //   Vector4 dq = rbd_.quaternionDivision(qtot1, qtot2);
  //   Vector3 da = rbd_.quaternion2mrp(dq);
  //   return da;
  // }

  // compute	the	control	input	using	w_t	=	x_{t+1}	-
  // \int_tˆ{t+1}f(x_t)
  Vector computeStateChange(DynPose& prev, double dt, InertiaRatios& ir) {
    Vector w;

    DynPose predicted = prev.propagate(dt, ir);

    Vector4 qTot = this->qTotal();
    Vector4 qTotpred = predicted.qTotal();
    Vector3 da = getMRPdifference(qTot, qTotpred);

    w = this->x() - predicted.x();
    w.segment<3>(6) = da;

    return w;
  }

  Vector6 getOdometry() {
    Vector6 odo;
    Vector x = rbd_.x();
    Vector4 qref = rbd_.qref();
    Vector3 a = x.segment<3>(6);
    odo.head(3) = x.segment<3>(0);
    Vector4 qnew = rbd_.addQuaternionError(a, qref);
    odo.tail(3) = rbd_.quaternion2mrp(qnew);
    return odo;
  }

  Vector6 getOdometry(DynPose& prev) {
    Vector6 odo;
    Vector x, xprev;
    Vector4 q, qprev;
    Vector3 a, aprev;

    x = rbd_.x();
    xprev = prev.x();
    a = x.segment<3>(6);

    q = rbd_.qref();

    qprev = prev.q();

    aprev = xprev.segment<3>(6);

    Vector3 dr = x.segment<3>(0) - xprev.segment<3>(0);
    Vector4 qtot_this = rbd_.addQuaternionError(a, q);
    Vector4 qtot_prev = rbd_.addQuaternionError(aprev, qprev);

    Vector4 qprev_inv;
    qprev_inv << -qtot_prev(0), -qtot_prev(1), -qtot_prev(2), qtot_prev(3);
    Vector4 qDiff = rbd_.quaternionMultiplication(qtot_this, qprev_inv);
    Vector3 mrp = rbd_.quaternion2mrp(qDiff);
    odo.tail(3) = mrp;

    Matrix3 Rprev = rbd_.rotationMatrix(qtot_prev);
    odo.head(3) = Rprev.transpose() * dr;

    return odo;
  }

  /* ************************************************************************ */

  DynPose getOdometryPose(DynPose& prev, bool initVelocities, double dt) {
    DynPose newPose(prev.rbd_.getIR(), prev.rbd_.getSigmaV(),
                    prev.rbd_.getSigmaW());
    Vector new_x(12);
    Vector3 new_r;
    Vector3 new_v;
    Vector3 new_a;
    Vector4 new_q;
    Vector3 new_w;

    Vector x, xprev;
    Vector4 q, qprev;
    Vector3 a, aprev;

    // get	x's
    x = rbd_.x();
    xprev = prev.x();

    // attitude	gets
    a = x.segment<3>(6);
    aprev = xprev.segment<3>(6);
    q = rbd_.qref();
    qprev = prev.q();
    // total	attitude
    Vector4 qtot_this = rbd_.addQuaternionError(a, q);
    Vector4 qtot_prev = rbd_.addQuaternionError(aprev, qprev);
    Vector4 qprev_inv;
    qprev_inv << -qtot_prev(0), -qtot_prev(1), -qtot_prev(2), qtot_prev(3);
    Vector4 qDiff = rbd_.quaternionMultiplication(qtot_this, qprev_inv);
    // previous	rotation	mat
    Matrix3 Rprev = rbd_.rotationMatrix(qtot_prev);

    new_r = Rprev.transpose() * (x.segment<3>(0) - xprev.segment<3>(0));
    Matrix3 Rdiff = rbd_.rotationMatrix(qDiff);
    new_q = rbd_.quaternionFromRot(Rdiff);

    if (initVelocities) {
      // differentiate	linear	velocity
      new_v = new_r / dt;

      /*	Differentiate	quaternion:
       *	dq/dt	=	(q[k]	-	q[k-1])/dt	= 0.5
       *O(w[k-1])	q[k-1]
       *	where	O(w[k-1])	is	an	orthonormal
       *quaternion	mult	matrix	for	[w1;	w2;
       *w3;	0]	(i.e.	quaternionMultiplication()) set	q[k-1]	=
       *[0;0;0;1]	by	definition	(from	a refererence frame)
       *and	solve	for	w[k-1]	gives
       *	w[k-1]	=	2	[q1[k];	q2[k];	q3[k]]	/	dt
       */
      new_w = 2 * new_q.head(3) / dt;
    } else {
      new_v = Vector3::Zero();
      new_w = Vector3::Zero();
    }
    new_a = Vector3::Zero();

    new_x.block<3, 1>(0, 0) = new_r;
    new_x.block<3, 1>(3, 0) = new_v;
    new_x.block<3, 1>(6, 0) = new_a;
    new_x.block<3, 1>(9, 0) = new_w;
    newPose.rbd_.setState(new_x, new_q);
    return newPose;
  }

  /* ************************************************************************ */

  DynPose adjustAttitude(DynPose& prev) {
    Vector4 q, qprev;
    DynPose newPose(prev.rbd_.getIR(), prev.rbd_.getSigmaV(),
                    prev.rbd_.getSigmaW());

    Vector x = rbd_.x();
    q = rbd_.qTotal();
    qprev = prev.qTotal();

    std::cout << "q:	" << q.transpose() << std::endl;
    std::cout << "qprev:	" << qprev.transpose() << std::endl;

    Matrix3 R = rbd_.rotationMatrix(q);
    Matrix3 Rprev = rbd_.rotationMatrix(qprev);
    Matrix3 Rdiff = R * Rprev.transpose();
    Vector4 new_qdiff = rbd_.quaternionFromRot(Rdiff);

    std::cout << "R:	" << R << std::endl;
    std::cout << "Rprev:	" << Rprev << std::endl;
    std::cout << "Rdiff:	" << Rdiff << std::endl;
    std::cout << "new_qdiff:	" << new_qdiff.transpose() << std::endl;

    Vector4 qnew = rbd_.quaternionMultiplication(new_qdiff, qprev);

    std::cout << "qnew	aa:	" << qnew.transpose() << std::endl << std::endl;
    if (isnan(qnew(1))) {
      std::cout << "qnew	aa	nan\n";
      new_qdiff = rbd_.quaternionFromRot(Rdiff);
    }

    x.segment<3>(6) = Vector3::Zero();
    rbd_.setState(x, qnew);
    newPose.rbd_.setState(x, qnew);
    return newPose;
  }

  /* ************************************************************************ */

  void shortenQuaternion(DynPose& prev) {
    Vector4 q, qprev, qnew;

    Vector x = rbd_.x();
    q = rbd_.qTotal();
    qprev = prev.qTotal();
    if (q.dot(qprev) < 0) {
      qnew = -q;
      x.segment<3>(6) = Vector3::Zero();
      rbd_.setState(x, qnew);
    }
  }

  /* ************************************************************************ */

  DynPose applyOdometry(DynPose& prev) {
    DynPose newPose(prev.rbd_.getIR(), prev.rbd_.getSigmaV(),
                    prev.rbd_.getSigmaW());
    Vector new_x(12);
    Vector3 new_r;
    Vector3 new_v;
    Vector3 new_a;
    Vector4 new_q;
    Vector3 new_w;

    Vector x, xprev;
    Vector4 q, qprev;
    Vector3 a, aprev;

    // get	x's
    x = rbd_.x();
    xprev = prev.x();

    // attitude	gets
    q = rbd_.qTotal();
    qprev = prev.qTotal();

    new_q = rbd_.quaternionMultiplication(q, qprev);

    Matrix3 Rprev = rbd_.rotationMatrix(qprev);
    new_r = Rprev * x.head(3) + xprev.head(3);

    new_v = Vector3::Zero();
    new_a = Vector3::Zero();
    new_w = Vector3::Zero();

    new_x.block<3, 1>(0, 0) = new_r;
    new_x.block<3, 1>(3, 0) = new_v;
    new_x.block<3, 1>(6, 0) = new_a;
    new_x.block<3, 1>(9, 0) = new_w;

    newPose.rbd_.setState(new_x, new_q);
    return newPose;
  }

  /* *********************************************************************** */

  Matrix4 wTo() const {
    Matrix4 T;

    // error	quaternion	is	applied
    Vector4 qtot = rbd_.qTotal();
    Vector x = rbd_.x();
    T.topLeftCorner(3, 3) = rbd_.rotationMatrix(qtot).transpose();
    T.col(3).head(3) << x.segment<3>(0);
    T.row(3) << 0., 0., 0., 1.;
    return T;
  }

  /* *********************************************************************** */

  Matrix4 oTw() const {
    Matrix4 T;
    Matrix3 R;

    // error	quaternion	is	applied
    Vector4 qtot = rbd_.qTotal();
    Vector x = rbd_.x();
    R = rbd_.rotationMatrix(qtot);

    T.topLeftCorner(3, 3) = R;
    T.col(3).head(3) << -R * x.segment<3>(0);
    T.row(3) << 0., 0., 0., 1.;
    return T;
  }

  Pose3 getPose3d() {
    return Pose3(
        this->wTo());  // may	be	wrong:	Mar	25,	2013,	B.E.T.
    // return	Pose3d(this->oTw());
  }

  /* *********************************************************************** */

  Vector4 transform_to_inertial(const Vector4& pBody) const {
    Vector3 p;
    p << pBody.x(), pBody.y(), pBody.z();
    Vector4 qtot = rbd_.qTotal();
    Vector x = rbd_.x();
    Vector3 T = x.head(3);
    Matrix3 Rt = rbd_.rotationMatrix(qtot).transpose();

    Vector3 pInertial = Rt * p + T;

    return Vector4(pInertial(0), pInertial(1), pInertial(2), 1.0);
  }

  /* *********************************************************************** */

  Vector4 transform_to_body(const Vector4& pInertial) const {
    Vector3 p;
    p << pInertial.x(), pInertial.y(), pInertial.z();
    Vector4 qtot = rbd_.qTotal();
    Vector x = rbd_.x();
    Vector3 T = x.head(3);
    Matrix3 R = rbd_.rotationMatrix(qtot);

    Vector3 pBody = R * (p - T);

    return Vector4(pBody(0), pBody(1), pBody(2), 1.0);
  }

  /* *********************************************************************** */

  Noise getProcessNoise(double dt, InertiaRatios ir) {
    Vector x0 = Vector::Zero(90);
    x0.head(12) = rbd_.x();
    rbd_.setIR(ir);
    Vector newLambda = rbd_.propagateRK4_adaptive(dt, x0).tail(78);

    Eigen::Matrix<double, 12, 12> lambda = rbd_.vec2symmMat(newLambda);
    Noise n = Covariance(lambda);
    return n;
  }

  /* *********************************************************************** */

  Vector vectorFull() const {
    Vector x = rbd_.x();
    Vector4 q = rbd_.qref();
    Vector3 mrp = rbd_.quaternion2mrp(q);
    x(6) += mrp(0);
    x(7) += mrp(1);
    x(8) += mrp(2);
    return x;
  }

  // Vector vector() const { return rbd_.x(); }

  RigidBody GetRBD() const { return rbd_; }

  /* ********************************************************************** */
  void print(const string& s) const { rbd_.print("Dynamic Pose:\n"); }

  void write(std::ostream& out) const {
    out << rbd_.print("Dynamic Pose:	\n ") << std::endl;
  }

  /* *********************************************************************** */
  bool equals(const DynPose& dp, double tol = 1e-9) const {
    return rbd_.equals(dp.GetRBD(), tol);
  }
};

std::ostream& operator<<(std::ostream& out, const DynPose& p) {
  p.write(out);
  return out;
}

}  // namespace gtsam
