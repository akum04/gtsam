#pragma once

// #include <isam/Factor.h>
// #include <isam/Node.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>

#include <gtsam/base/Testable.h>
#include <gtsam/nonlinear/NonlinearFactor.h>

#include <string>

// #include "NodeExmap.h"

namespace gtsam {

class InertiaRatios {
  double _k1, _k2;

  friend std::ostream& operator<<(std::ostream& out, const InertiaRatios& p) {
    p.write(out);
    return out;
  }

  /*
   *	k1	=	ln(J11	/	J22)
   *	k2	=	ln(J22	/	J33)
   */

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  static const int dim = 2;
  static const char* name() { return "InertiaRatios"; }

  InertiaRatios() {
    _k1 = 0;
    _k2 = 0;
  }

  InertiaRatios(Vector2& k) {
    _k1 = k(0);
    _k2 = k(1);
  }

  InertiaRatios(const double& k1, const double& k2) {
    _k1 = k1;
    _k2 = k2;
  }

  Matrix3 getJ() const {
    Matrix3 J = Matrix3::Zero();
    double Jscale = 1.0;  // 0.0116;
    J(0, 0) = exp(_k1);
    J(1, 1) = 1.0;
    J(2, 2) = exp(-_k2);

    J *= Jscale;

    return J;
  }

  Eigen::VectorXd x() const {
    Vector2 x;
    x(0) = _k1;
    x(1) = _k2;
    return x;
  }

  void setState(Eigen::VectorXd x) {
    _k1 = x(0);
    _k2 = x(1);
  }

  InertiaRatios exmap(const Vector2& Delta) {
    InertiaRatios res = *this;
    res._k1 += Delta(0);
    res._k2 += Delta(1);
    return res;
  }

  InertiaRatios exmap_reset(const Vector2& Delta) {
    InertiaRatios res = *this;
    res._k1 += Delta(0);
    res._k2 += Delta(1);
    return res;
  }

  Eigen::VectorXd vector() const {
    Vector2 tmp;
    tmp << _k1, _k2;
    return tmp;
  }

  void set(const Vector2& v) {
    _k1 = v(0);
    _k2 = v(1);
  }

  void write(std::ostream& out) const {
    Matrix3 Jcurr = getJ();
    out << std::endl
        << "inertaRatios	x:	" << x().transpose() << std::endl
        << Jcurr(0, 0) << "	,	" << Jcurr(1, 1) << "	,	"
        << Jcurr(2, 2) << std::endl;
  }

  bool equals(const InertiaRatios& ir, double tol = 1e-9) const {
    return (_k1 < tol && _k2 < tol);
  }
};

// typedef NodeExmapT<InertiaRatios> InertiaRatios_Node;

// class InertiaRatiosFactor : gtsam::NoiseModelFactor1<InertiaRatios> {
//  public:
//   Eigen::VectorXd evaluateError(Selector s = ESTIMATE) const {
//     InertiaRatios ir = _ir_node->value(s);
//     Eigen::VectorXd err = ir.vector() - _measure.vector();
//     return err;
//   }
// };
/**
 *	Prior	on	InertiaRatios.
 */
/*
class InertiaRatios_Factor : public FactorT<InertiaRatios> {
 public:
  InertiaRatios_Node* _ir_node;

  InertiaRatios_Factor(InertiaRatios_Node* ir_node, const InertiaRatios& prior,
                       const Noise& noise)
      : FactorT<InertiaRatios>("InertiaRatios_Factor", 2, noise, prior),
        _ir_node(ir_node) {
    _nodes.resize(1);
    _nodes[0] = ir_node;
  }

  void initialize() {
    if (!_ir_node->initialized()) {
      InertiaRatios predict = _measure;
      _ir_node->init(predict);
    }
  }

  Eigen::VectorXd basic_error(Selector s = ESTIMATE) const {
    InertiaRatios ir = _ir_node->value(s);
    Eigen::VectorXd err = ir.vector() - _measure.vector();
    return err;
  }
};
*/
}  // namespace gtsam
