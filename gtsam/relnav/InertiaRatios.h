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
  double k1_, k2_;

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
    k1_ = 0;
    k2_ = 0;
  }

  InertiaRatios(Vector2& k) {
    k1_ = k(0);
    k2_ = k(1);
  }

  InertiaRatios(const double& k1, const double& k2) {
    k1_ = k1;
    k2_ = k2;
  }

  Matrix3 getJ() const {
    Matrix3 J = Matrix3::Zero();
    double Jscale = 1.0;  // 0.0116;
    J(0, 0) = exp(k1_);
    J(1, 1) = 1.0;
    J(2, 2) = exp(-k2_);

    J *= Jscale;

    return J;
  }

  Vector2 x() const {
    Vector2 x;
    x(0) = k1_;
    x(1) = k2_;
    return x;
  }

  void setState(Vector2 x) {
    k1_ = x(0);
    k2_ = x(1);
  }

  InertiaRatios exmap(const Vector2& Delta) {
    InertiaRatios res = *this;
    res.k1_ += Delta(0);
    res.k2_ += Delta(1);
    return res;
  }

  InertiaRatios exmap_reset(const Vector2& Delta) {
    InertiaRatios res = *this;
    res.k1_ += Delta(0);
    res.k2_ += Delta(1);
    return res;
  }

  Vector2 vector() const {
    Vector2 tmp;
    tmp << k1_, k2_;
    return tmp;
  }

  void set(const Vector2& v) {
    k1_ = v(0);
    k2_ = v(1);
  }

  void write(std::ostream& out) const {
    Matrix3 Jcurr = getJ();
    out << std::endl
        << "inertaRatios	x:	" << x().transpose() << std::endl
        << Jcurr(0, 0) << "	,	" << Jcurr(1, 1) << "	,	"
        << Jcurr(2, 2) << std::endl;
  }

  bool equals(const InertiaRatios& ir, double tol = 1e-9) const {
    return (k1_ < tol && k2_ < tol);
  }

  void print(const std::string& str) {
    std::cout << k1_ << "\n" << k2_ << std::endl;
  }
};
/*
typedef NodeExmapT<InertiaRatios> InertiaRatios_Node;

class InertiaRatiosFactor : gtsam::NoiseModelFactor1<InertiaRatios> {
 public:
  Vector evaluateError(Selector s = ESTIMATE) const {
    InertiaRatios ir = _ir_node->value(s);
    Vector err = ir.vector() - _measure.vector();
    return err;
  }
};
*/

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

  Vector basic_error(Selector s = ESTIMATE) const {
    InertiaRatios ir = _ir_node->value(s);
    Vector err = ir.vector() - _measure.vector();
    return err;
  }
};
*/
}  // namespace gtsam
