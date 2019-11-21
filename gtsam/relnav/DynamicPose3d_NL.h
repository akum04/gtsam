#pragma once

#include <Eigen/Dense>
#include <ostream>
// #include "FactorVariableNoise.h"
// #include "isam/isam.h"
#include <gtsam/geometry/Pose3.h>
#include "Noise.h"
#include "RigidBodyDynamics.h"

using namespace Eigen;
namespace gtsam {

class dynamicPose3d_NL {
  friend std::ostream& operator<<(std::ostream& out,
                                  const dynamicPose3d_NL& p) {
    p.write(out);
    return out;
  }

  RigidBodyDynamics rbd;

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  //	assignment	operator	and	copy	constructor
  // implicitly created,	which	is	ok
  static const int dim = 12;
  static const char* name() { return "dynamicPose3d_NL"; }

  Noise* factor_noise;  // check	if	this	is	ever	used

  dynamicPose3d_NL(InertiaRatios ir, double sigma_v, double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {}

  // copy	constructor
  dynamicPose3d_NL(const dynamicPose3d_NL& cSource)
      : rbd(cSource.rbd.getIR(), cSource.rbd.getSigmaV(),
            cSource.rbd.getSigmaW()) {
    rbd.setState(cSource.rbd.x(), cSource.rbd.qref());
  }

  dynamicPose3d_NL& operator=(const dynamicPose3d_NL& cSource) {
    rbd = RigidBodyDynamics(cSource.rbd.getIR(), cSource.rbd.getSigmaV(),
                            cSource.rbd.getSigmaW());
    rbd.setState(cSource.rbd.x(), cSource.rbd.qref());
    return *this;
  }

  dynamicPose3d_NL(VectorXd x, InertiaRatios ir, double sigma_v, double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {
    Vector4d qref;
    qref << 0, 0, 0, 1;
    if (x.size() == 12) {
      rbd.setState(x, qref);
    }
  }

  dynamicPose3d_NL(VectorXd x, Vector4d qref, InertiaRatios ir, double sigma_v,
                   double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {
    if (x.size() == 12) {
      rbd.setState(x, qref);
    }
  }

  dynamicPose3d_NL(const Matrix4d& hm, bool initVelocities, double dt,
                   InertiaRatios ir, double sigma_v, double sigma_w)
      : rbd(ir, sigma_v, sigma_w) {
    Eigen::Matrix<double, 12, 1> x;
    Vector3d r;
    Vector3d v;
    Vector3d a;
    Vector4d q;
    Vector3d w;

    // Convert	matrix	to	R,T
    Matrix4d HM = hm / hm(3, 3);  //	enforce	T(3,3)=1
    Matrix3d R = HM.topLeftCorner(3, 3);
    Vector3d T = HM.col(3).head(3);

    // compute	quaternion
    q = rbd.quaternionFromRot(R);
    a = Vector3d::Zero();

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
      v = Vector3d::Zero();
      w = Vector3d::Zero();
    }

    x.block<3, 1>(0, 0) = T;
    x.block<3, 1>(3, 0) = v;
    x.block<3, 1>(6, 0) = a;
    x.block<3, 1>(9, 0) = w;
    rbd.setState(x, q);
  }

  VectorXd x() { return rbd.x(); };
  Vector4d q() { return rbd.qref(); };
  Vector4d qTotal() const { return rbd.qTotal(); };

  dynamicPose3d_NL exmap(const Eigen::Matrix<double, 12, 1>& Delta) const {
    dynamicPose3d_NL res = *this;
    res.rbd.setState(res.rbd.x() + Delta);
    return res;
  }

  dynamicPose3d_NL exmap_reset(const Eigen::Matrix<double, 12, 1>& Delta) {
    dynamicPose3d_NL res = *this;
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

  void set(const VectorXd& v) { rbd.setState(v); }

  void set_qref(const Vector4d& qset) { rbd.setState(rbd.x(), qset); }

  void rezero() {
    VectorXd x = VectorXd::Zero(12);
    Vector4d q;
    q << 0, 0, 0, 1;
    rbd.setState(x, q);
  }

  dynamicPose3d_NL propagate(double dt, InertiaRatios& ir) {
    VectorXd x0 = VectorXd::Zero(90);
    x0.head(12) = rbd.x();
    rbd.setIR(ir);
    //	std::cout	<<	"x0:	"	<<	x0.head(12).transpose()
    //<< std::endl;
    VectorXd newX = rbd.propagateRK4_adaptive(dt, x0).head(12);

    //	std::cout	<<	"dt:	"	<<	dt	<<
    // std::endl; 	std::cout	<<	"newX:	"	<<
    // newX.transpose()
    //<<	std::endl;

    dynamicPose3d_NL xNew(newX, this->rbd.qref(), this->rbd.getIR(),
                          this->rbd.getSigmaV(), this->rbd.getSigmaW());
    xNew.exmap(Eigen::Matrix<double, 12, 1>::Zero());
    return xNew;
  }

  Vector3d getMRPdifference(Vector4d qtot1, Vector4d qtot2) {
    Vector4d dq = rbd.quaternionDivision(qtot1, qtot2);
    Vector3d da = rbd.quaternion2mrp(dq);
    return da;
  }

  // compute	the	control	input	using	w_t	=	x_{t+1}	-
  // \int_tˆ{t+1}f(x_t)
  VectorXd computeStateChange(dynamicPose3d_NL& prev, double dt,
                              InertiaRatios& ir) {
    VectorXd w;

    dynamicPose3d_NL predicted = prev.propagate(dt, ir);

    Vector4d qTot = this->qTotal();
    Vector4d qTotpred = predicted.qTotal();
    Vector3d da = getMRPdifference(qTot, qTotpred);

    w = this->x() - predicted.x();
    w.segment<3>(6) = da;

    return w;
  }

  gtsam::Vector6 getOdometry() {
    gtsam::Vector6 odo;
    VectorXd x = rbd.x();
    Vector4d qref = rbd.qref();
    Vector3d a = x.segment<3>(6);
    odo.head(3) = x.segment<3>(0);
    Vector4d qnew = rbd.addQuaternionError(a, qref);
    odo.tail(3) = rbd.quaternion2mrp(qnew);
    return odo;
  }

  gtsam::Vector6 getOdometry(dynamicPose3d_NL& prev) {
    gtsam::Vector6 odo;
    VectorXd x, xprev;
    Vector4d q, qprev;
    Vector3d a, aprev;

    x = rbd.x();
    xprev = prev.x();
    a = x.segment<3>(6);

    q = rbd.qref();

    qprev = prev.q();

    aprev = xprev.segment<3>(6);

    Vector3d dr = x.segment<3>(0) - xprev.segment<3>(0);
    Vector4d qtot_this = rbd.addQuaternionError(a, q);
    Vector4d qtot_prev = rbd.addQuaternionError(aprev, qprev);

    Vector4d qprev_inv;
    qprev_inv << -qtot_prev(0), -qtot_prev(1), -qtot_prev(2), qtot_prev(3);
    Vector4d qDiff = rbd.quaternionMultiplication(qtot_this, qprev_inv);
    Vector3d mrp = rbd.quaternion2mrp(qDiff);
    odo.tail(3) = mrp;

    Matrix3d Rprev = rbd.rotationMatrix(qtot_prev);
    odo.head(3) = Rprev.transpose() * dr;

    return odo;
  }

  dynamicPose3d_NL getOdometryPose(dynamicPose3d_NL& prev, bool initVelocities,
                                   double dt) {
    dynamicPose3d_NL newPose(prev.rbd.getIR(), prev.rbd.getSigmaV(),
                             prev.rbd.getSigmaW());
    VectorXd new_x(12);
    Vector3d new_r;
    Vector3d new_v;
    Vector3d new_a;
    Vector4d new_q;
    Vector3d new_w;

    VectorXd x, xprev;
    Vector4d q, qprev;
    Vector3d a, aprev;

    // get	x's
    x = rbd.x();
    xprev = prev.x();

    // attitude	gets
    a = x.segment<3>(6);
    aprev = xprev.segment<3>(6);
    q = rbd.qref();
    qprev = prev.q();
    // total	attitude
    Vector4d qtot_this = rbd.addQuaternionError(a, q);
    Vector4d qtot_prev = rbd.addQuaternionError(aprev, qprev);
    Vector4d qprev_inv;
    qprev_inv << -qtot_prev(0), -qtot_prev(1), -qtot_prev(2), qtot_prev(3);
    Vector4d qDiff = rbd.quaternionMultiplication(qtot_this, qprev_inv);
    // previous	rotation	mat
    Matrix3d Rprev = rbd.rotationMatrix(qtot_prev);

    new_r = Rprev.transpose() * (x.segment<3>(0) - xprev.segment<3>(0));
    Matrix3d Rdiff = rbd.rotationMatrix(qDiff);
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
      new_v = Vector3d::Zero();
      new_w = Vector3d::Zero();
    }
    new_a = Vector3d::Zero();

    new_x.block<3, 1>(0, 0) = new_r;
    new_x.block<3, 1>(3, 0) = new_v;
    new_x.block<3, 1>(6, 0) = new_a;
    new_x.block<3, 1>(9, 0) = new_w;
    newPose.rbd.setState(new_x, new_q);
    return newPose;
  }

  dynamicPose3d_NL adjustAttitude(dynamicPose3d_NL& prev) {
    Vector4d q, qprev;
    dynamicPose3d_NL newPose(prev.rbd.getIR(), prev.rbd.getSigmaV(),
                             prev.rbd.getSigmaW());

    VectorXd x = rbd.x();
    q = rbd.qTotal();
    qprev = prev.qTotal();

    std::cout << "q:	" << q.transpose() << std::endl;
    std::cout << "qprev:	" << qprev.transpose() << std::endl;

    Matrix3d R = rbd.rotationMatrix(q);
    Matrix3d Rprev = rbd.rotationMatrix(qprev);
    Matrix3d Rdiff = R * Rprev.transpose();
    Vector4d new_qdiff = rbd.quaternionFromRot(Rdiff);

    std::cout << "R:	" << R << std::endl;
    std::cout << "Rprev:	" << Rprev << std::endl;
    std::cout << "Rdiff:	" << Rdiff << std::endl;
    std::cout << "new_qdiff:	" << new_qdiff.transpose() << std::endl;

    Vector4d qnew = rbd.quaternionMultiplication(new_qdiff, qprev);

    std::cout << "qnew	aa:	" << qnew.transpose() << std::endl << std::endl;
    if (isnan(qnew(1))) {
      std::cout << "qnew	aa	nan\n";
      new_qdiff = rbd.quaternionFromRot(Rdiff);
    }

    x.segment<3>(6) = Vector3d::Zero();
    rbd.setState(x, qnew);
    newPose.rbd.setState(x, qnew);
    return newPose;
  }

  void shortenQuaternion(dynamicPose3d_NL& prev) {
    Vector4d q, qprev, qnew;

    VectorXd x = rbd.x();
    q = rbd.qTotal();
    qprev = prev.qTotal();
    if (q.dot(qprev) < 0) {
      qnew = -q;
      x.segment<3>(6) = Vector3d::Zero();
      rbd.setState(x, qnew);
    }
  }

  dynamicPose3d_NL applyOdometry(dynamicPose3d_NL& prev) {
    dynamicPose3d_NL newPose(prev.rbd.getIR(), prev.rbd.getSigmaV(),
                             prev.rbd.getSigmaW());
    VectorXd new_x(12);
    Vector3d new_r;
    Vector3d new_v;
    Vector3d new_a;
    Vector4d new_q;
    Vector3d new_w;

    VectorXd x, xprev;
    Vector4d q, qprev;
    Vector3d a, aprev;

    // get	x's
    x = rbd.x();
    xprev = prev.x();

    // attitude	gets
    q = rbd.qTotal();
    qprev = prev.qTotal();

    new_q = rbd.quaternionMultiplication(q, qprev);

    Matrix3d Rprev = rbd.rotationMatrix(qprev);
    new_r = Rprev * x.head(3) + xprev.head(3);

    new_v = Vector3d::Zero();
    new_a = Vector3d::Zero();
    new_w = Vector3d::Zero();

    new_x.block<3, 1>(0, 0) = new_r;
    new_x.block<3, 1>(3, 0) = new_v;
    new_x.block<3, 1>(6, 0) = new_a;
    new_x.block<3, 1>(9, 0) = new_w;

    newPose.rbd.setState(new_x, new_q);
    return newPose;
  }

  Matrix4d wTo() const {
    Matrix4d T;

    // error	quaternion	is	applied
    Vector4d qtot = rbd.qTotal();
    VectorXd x = rbd.x();
    T.topLeftCorner(3, 3) = rbd.rotationMatrix(qtot).transpose();
    T.col(3).head(3) << x.segment<3>(0);
    T.row(3) << 0., 0., 0., 1.;
    return T;
  }

  Matrix4d oTw() const {
    Matrix4d T;
    Matrix3d R;

    // error	quaternion	is	applied
    Vector4d qtot = rbd.qTotal();
    VectorXd x = rbd.x();
    R = rbd.rotationMatrix(qtot);

    T.topLeftCorner(3, 3) = R;
    T.col(3).head(3) << -R * x.segment<3>(0);
    T.row(3) << 0., 0., 0., 1.;
    return T;
  }

  gtsam::Pose3 getPose3d() {
    return gtsam::Pose3(
        this->wTo());  // may	be	wrong:	Mar	25,	2013,	B.E.T.
    // return	Pose3d(this->oTw());
  }

  gtsam::Vector4 transform_to_inertial(const gtsam::Vector4& pBody) const {
    Vector3d p;
    p << pBody.x(), pBody.y(), pBody.z();
    Vector4d qtot = rbd.qTotal();
    VectorXd x = rbd.x();
    Vector3d T = x.head(3);
    Matrix3d Rt = rbd.rotationMatrix(qtot).transpose();

    Vector3d pInertial = Rt * p + T;

    return gtsam::Vector4(pInertial(0), pInertial(1), pInertial(2), 1.0);
  }

  gtsam::Vector4 transform_to_body(const gtsam::Vector4& pInertial) const {
    Vector3d p;
    p << pInertial.x(), pInertial.y(), pInertial.z();
    Vector4d qtot = rbd.qTotal();
    VectorXd x = rbd.x();
    Vector3d T = x.head(3);
    Matrix3d R = rbd.rotationMatrix(qtot);

    Vector3d pBody = R * (p - T);

    return gtsam::Vector4(pBody(0), pBody(1), pBody(2), 1.0);
  }

  gtsam::Noise getProcessNoise(double dt, InertiaRatios ir) {
    VectorXd x0 = VectorXd::Zero(90);
    x0.head(12) = rbd.x();
    rbd.setIR(ir);
    VectorXd newLambda = rbd.propagateRK4_adaptive(dt, x0).tail(78);

    Eigen::Matrix<double, 12, 12> lambda = rbd.vec2symmMat(newLambda);
    gtsam::Noise n = gtsam::Covariance(lambda);
    return n;
  }

  VectorXd vectorFull() const {
    VectorXd x = rbd.x();
    Vector4d q = rbd.qref();
    Vector3d mrp = rbd.quaternion2mrp(q);
    x(6) += mrp(0);
    x(7) += mrp(1);
    x(8) += mrp(2);
    return x;
  }

  VectorXd vector() const { return rbd.x(); }

  void write(std::ostream& out) const {
    out << std::endl
        << "dP3d_NL	x:	" << rbd.x().transpose() << std::endl;
    out << "dP3d_NL	qref:	" << rbd.qref().transpose() << std::endl;
    out << std::endl;
  }
};
}  // namespace gtsam
