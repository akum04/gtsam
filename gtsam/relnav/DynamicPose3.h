#pragma once

#include <Eigen/Dense>
#include <ostream>
// #include "FactorVariableNoise.h"
// #include "isam/isam.h"
#include <gtsam/geometry/Pose3.h>
#include "Noise.h"
#include "RigidBodyDynamics.h"

// using namespace Eigen;
namespace gtsam {

class DynamicPose3 {
 private:
  RigidBodyDynamics rbd;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  //	assignment	operator	and	copy	constructor
  // implicitly created,	which	is	ok
  static const int dim = 12;
  static const char* name() { return "DynamicPose3"; }

  Noise* factor_noise;  // check	if	this	is	ever	used

  DynamicPose3(InertiaRatios ir, double sigma_v, double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {}

  // copy	constructor
  DynamicPose3(const DynamicPose3& cSource)
      : rbd(cSource.rbd.getIR(), cSource.rbd.getSigmaV(),
            cSource.rbd.getSigmaW()) {
    rbd.setState(cSource.rbd.x(), cSource.rbd.qref());
  }

  DynamicPose3& operator=(const DynamicPose3& cSource) {
    rbd = RigidBodyDynamics(cSource.rbd.getIR(), cSource.rbd.getSigmaV(),
                            cSource.rbd.getSigmaW());
    rbd.setState(cSource.rbd.x(), cSource.rbd.qref());
    return *this;
  }

  DynamicPose3(Vector x, InertiaRatios ir, double sigma_v, double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {
    Vector4 qref;
    qref << 0, 0, 0, 1;
    if (x.size() == 12) {
      rbd.setState(x, qref);
    }
  }

  DynamicPose3(Vector x, Vector4 qref, InertiaRatios ir, double sigma_v,
               double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {
    if (x.size() == 12) {
      rbd.setState(x, qref);
    }
  }

  DynamicPose3(const Matrix4& hm, bool initVelocities, double dt,
               InertiaRatios ir, double sigma_v, double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {
    Vector12 x;
    Vector3 r;
    Vector3 v;
    Vector3 a;
    Vector4 q;
    Vector3 w;

    // Convert	matrix	to	R,T
    Matrix4 HM = hm / hm(3, 3);  //	enforce	T(3,3)=1
    Matrix3 R = HM.topLeftCorner(3, 3);
    Vector3 T = HM.col(3).head(3);

    // compute	quaternion
    q = rbd.quaternionFromRot(R);
    a = Vector3::Zero();

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
      w = 2 * q.head(3) / dt;
    } else {
      v = Vector3::Zero();
      w = Vector3::Zero();
    }

    x.block<3, 1>(0, 0) = T;
    x.block<3, 1>(3, 0) = v;
    x.block<3, 1>(6, 0) = a;
    x.block<3, 1>(9, 0) = w;
    rbd.setState(x, q);
  }

  Vector x() { return rbd.x(); };
  Vector4 q() { return rbd.qref(); };
  Vector4 qTotal() const { return rbd.qTotal(); };

  DynamicPose3 exmap(const Vector12& Delta) const {
    DynamicPose3 res = *this;
    res.rbd.setState(res.rbd.x() + Delta);
    return res;
  }

  DynamicPose3 exmap_reset(const Vector12& Delta) {
    DynamicPose3 res = *this;
    res.rbd.setState(res.rbd.x() + Delta);
    res.rbd.reset_qref();

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
    // NodeExmapT->apply_reset();

    return res;
  }

  void set(const Vector& v) { rbd.setState(v); }

  void set_qref(const Vector4& qset) { rbd.setState(rbd.x(), qset); }

  void rezero() {
    Vector x = Vector::Zero(12);
    Vector4 q;
    q << 0, 0, 0, 1;
    rbd.setState(x, q);
  }

  DynamicPose3 propagate(double dt, InertiaRatios& ir) {
    Vector x0 = Vector::Zero(90);
    x0.head(12) = rbd.x();
    rbd.setIR(ir);
    //	std::cout	<<	"x0:	"	<<	x0.head(12).transpose()
    //<< std::endl;
    Vector newX = rbd.propagateRK4_adaptive(dt, x0).head(12);

    //	std::cout	<<	"dt:	"	<<	dt	<<
    // std::endl; 	std::cout	<<	"newX:	"	<<
    // newX.transpose()
    //<<	std::endl;

    DynamicPose3 xNew(newX, this->rbd.qref(), this->rbd.getIR(),
                      this->rbd.getSigmaV(), this->rbd.getSigmaW());
    xNew.exmap(Vector12::Zero());
    return xNew;
  }

  Vector3 getMRPdifference(Vector4 qtot1, Vector4 qtot2) {
    Vector4 dq = rbd.quaternionDivision(qtot1, qtot2);
    Vector3 da = rbd.quaternion2mrp(dq);
    return da;
  }

  // compute	the	control	input	using	w_t	=	x_{t+1}	-
  // \int_tˆ{t+1}f(x_t)
  Vector computeStateChange(DynamicPose3& prev, double dt, InertiaRatios& ir) {
    Vector w;

    DynamicPose3 predicted = prev.propagate(dt, ir);

    Vector4 qTot = this->qTotal();
    Vector4 qTotpred = predicted.qTotal();
    Vector3 da = getMRPdifference(qTot, qTotpred);

    w = this->x() - predicted.x();
    w.segment<3>(6) = da;

    return w;
  }

  Vector6 getOdometry() {
    Vector6 odo;
    Vector x = rbd.x();
    Vector4 qref = rbd.qref();
    Vector3 a = x.segment<3>(6);
    odo.head(3) = x.segment<3>(0);
    Vector4 qnew = rbd.addQuaternionError(a, qref);
    odo.tail(3) = rbd.quaternion2mrp(qnew);
    return odo;
  }

  Vector6 getOdometry(DynamicPose3& prev) {
    Vector6 odo;
    Vector x, xprev;
    Vector4 q, qprev;
    Vector3 a, aprev;

    x = rbd.x();
    xprev = prev.x();
    a = x.segment<3>(6);

    q = rbd.qref();

    qprev = prev.q();

    aprev = xprev.segment<3>(6);

    Vector3 dr = x.segment<3>(0) - xprev.segment<3>(0);
    Vector4 qtot_this = rbd.addQuaternionError(a, q);
    Vector4 qtot_prev = rbd.addQuaternionError(aprev, qprev);

    Vector4 qprev_inv;
    qprev_inv << -qtot_prev(0), -qtot_prev(1), -qtot_prev(2), qtot_prev(3);
    Vector4 qDiff = rbd.quaternionMultiplication(qtot_this, qprev_inv);
    Vector3 mrp = rbd.quaternion2mrp(qDiff);
    odo.tail(3) = mrp;

    Matrix3 Rprev = rbd.rotationMatrix(qtot_prev);
    odo.head(3) = Rprev.transpose() * dr;

    return odo;
  }

  DynamicPose3 getOdometryPose(DynamicPose3& prev, bool initVelocities,
                               double dt) {
    DynamicPose3 newPose(prev.rbd.getIR(), prev.rbd.getSigmaV(),
                         prev.rbd.getSigmaW());
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
    x = rbd.x();
    xprev = prev.x();

    // attitude	gets
    a = x.segment<3>(6);
    aprev = xprev.segment<3>(6);
    q = rbd.qref();
    qprev = prev.q();
    // total	attitude
    Vector4 qtot_this = rbd.addQuaternionError(a, q);
    Vector4 qtot_prev = rbd.addQuaternionError(aprev, qprev);
    Vector4 qprev_inv;
    qprev_inv << -qtot_prev(0), -qtot_prev(1), -qtot_prev(2), qtot_prev(3);
    Vector4 qDiff = rbd.quaternionMultiplication(qtot_this, qprev_inv);
    // previous	rotation	mat
    Matrix3 Rprev = rbd.rotationMatrix(qtot_prev);

    new_r = Rprev.transpose() * (x.segment<3>(0) - xprev.segment<3>(0));
    Matrix3 Rdiff = rbd.rotationMatrix(qDiff);
    new_q = rbd.quaternionFromRot(Rdiff);

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
    newPose.rbd.setState(new_x, new_q);
    return newPose;
  }

  DynamicPose3 adjustAttitude(DynamicPose3& prev) {
    Vector4 q, qprev;
    DynamicPose3 newPose(prev.rbd.getIR(), prev.rbd.getSigmaV(),
                         prev.rbd.getSigmaW());

    Vector x = rbd.x();
    q = rbd.qTotal();
    qprev = prev.qTotal();

    std::cout << "q:	" << q.transpose() << std::endl;
    std::cout << "qprev:	" << qprev.transpose() << std::endl;

    Matrix3 R = rbd.rotationMatrix(q);
    Matrix3 Rprev = rbd.rotationMatrix(qprev);
    Matrix3 Rdiff = R * Rprev.transpose();
    Vector4 new_qdiff = rbd.quaternionFromRot(Rdiff);

    std::cout << "R:	" << R << std::endl;
    std::cout << "Rprev:	" << Rprev << std::endl;
    std::cout << "Rdiff:	" << Rdiff << std::endl;
    std::cout << "new_qdiff:	" << new_qdiff.transpose() << std::endl;

    Vector4 qnew = rbd.quaternionMultiplication(new_qdiff, qprev);

    std::cout << "qnew	aa:	" << qnew.transpose() << std::endl << std::endl;
    if (isnan(qnew(1))) {
      std::cout << "qnew	aa	nan\n";
      new_qdiff = rbd.quaternionFromRot(Rdiff);
    }

    x.segment<3>(6) = Vector3::Zero();
    rbd.setState(x, qnew);
    newPose.rbd.setState(x, qnew);
    return newPose;
  }

  void shortenQuaternion(DynamicPose3& prev) {
    Vector4 q, qprev, qnew;

    Vector x = rbd.x();
    q = rbd.qTotal();
    qprev = prev.qTotal();
    if (q.dot(qprev) < 0) {
      qnew = -q;
      x.segment<3>(6) = Vector3::Zero();
      rbd.setState(x, qnew);
    }
  }

  DynamicPose3 applyOdometry(DynamicPose3& prev) {
    DynamicPose3 newPose(prev.rbd.getIR(), prev.rbd.getSigmaV(),
                         prev.rbd.getSigmaW());
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
    x = rbd.x();
    xprev = prev.x();

    // attitude	gets
    q = rbd.qTotal();
    qprev = prev.qTotal();

    new_q = rbd.quaternionMultiplication(q, qprev);

    Matrix3 Rprev = rbd.rotationMatrix(qprev);
    new_r = Rprev * x.head(3) + xprev.head(3);

    new_v = Vector3::Zero();
    new_a = Vector3::Zero();
    new_w = Vector3::Zero();

    new_x.block<3, 1>(0, 0) = new_r;
    new_x.block<3, 1>(3, 0) = new_v;
    new_x.block<3, 1>(6, 0) = new_a;
    new_x.block<3, 1>(9, 0) = new_w;

    newPose.rbd.setState(new_x, new_q);
    return newPose;
  }

  Matrix4 wTo() const {
    Matrix4 T;

    // error	quaternion	is	applied
    Vector4 qtot = rbd.qTotal();
    Vector x = rbd.x();
    T.topLeftCorner(3, 3) = rbd.rotationMatrix(qtot).transpose();
    T.col(3).head(3) << x.segment<3>(0);
    T.row(3) << 0., 0., 0., 1.;
    return T;
  }

  Matrix4 oTw() const {
    Matrix4 T;
    Matrix3 R;

    // error	quaternion	is	applied
    Vector4 qtot = rbd.qTotal();
    Vector x = rbd.x();
    R = rbd.rotationMatrix(qtot);

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

  Vector4 transform_to_inertial(const Vector4& pBody) const {
    Vector3 p;
    p << pBody.x(), pBody.y(), pBody.z();
    Vector4 qtot = rbd.qTotal();
    Vector x = rbd.x();
    Vector3 T = x.head(3);
    Matrix3 Rt = rbd.rotationMatrix(qtot).transpose();

    Vector3 pInertial = Rt * p + T;

    return Vector4(pInertial(0), pInertial(1), pInertial(2), 1.0);
  }

  Vector4 transform_to_body(const Vector4& pInertial) const {
    Vector3 p;
    p << pInertial.x(), pInertial.y(), pInertial.z();
    Vector4 qtot = rbd.qTotal();
    Vector x = rbd.x();
    Vector3 T = x.head(3);
    Matrix3 R = rbd.rotationMatrix(qtot);

    Vector3 pBody = R * (p - T);

    return Vector4(pBody(0), pBody(1), pBody(2), 1.0);
  }

  Noise getProcessNoise(double dt, InertiaRatios ir) {
    Vector x0 = Vector::Zero(90);
    x0.head(12) = rbd.x();
    rbd.setIR(ir);
    Vector newLambda = rbd.propagateRK4_adaptive(dt, x0).tail(78);

    Eigen::Matrix<double, 12, 12> lambda = rbd.vec2symmMat(newLambda);
    Noise n = Covariance(lambda);
    return n;
  }

  Vector vectorFull() const {
    Vector x = rbd.x();
    Vector4 q = rbd.qref();
    Vector3 mrp = rbd.quaternion2mrp(q);
    x(6) += mrp(0);
    x(7) += mrp(1);
    x(8) += mrp(2);
    return x;
  }

  Vector vector() const { return rbd.x(); }

  void write(std::ostream& out) const {
    out << std::endl
        << "dP3d_NL	x:	" << rbd.x().transpose() << std::endl;
    out << "dP3d_NL	qref:	" << rbd.qref().transpose() << std::endl;
    out << std::endl;
  }

  bool equals(const DynamicPose3& dp3, double tol = 1e-9) const { return true; }
};

std::ostream& operator<<(std::ostream& out, const DynamicPose3& p) {
  p.write(out);
  return out;
}

}  // namespace gtsam
