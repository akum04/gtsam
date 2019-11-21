/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  InertiaRatiosFactor.h
 *  @author Luca Carlone
 **/
#pragma once

#include <ostream>

#include <gtsam/base/VectorSpace.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include "InertiaRatios.h"

namespace gtsam {

/**
 * A class to model GPS measurements, including a bias term which models
 * common-mode errors and that can be partially corrected if other sensors are
 * used
 * @addtogroup SLAM
 */
// template <class InertiaRatios>
class InertiaRatiosFactor : public NoiseModelFactor1<InertiaRatios> {
 private:
  typedef InertiaRatiosFactor This;
  typedef NoiseModelFactor1<InertiaRatios> Base;

  InertiaRatios prior_; /** The measurement */

 public:
  // shorthand for a smart pointer to a factor
  typedef boost::shared_ptr<InertiaRatiosFactor> shared_ptr;

  /** default constructor - only use for serialization */
  InertiaRatiosFactor() {}

  /** Constructor */
  InertiaRatiosFactor(Key key, const InertiaRatios& prior,
                      const SharedNoiseModel& model)
      : Base(model, key), prior_(prior) {}

  virtual ~InertiaRatiosFactor() {}

  /** implement functions needed for Testable */

  /** print */
  virtual void print(const std::string& s, const KeyFormatter& keyFormatter =
                                               DefaultKeyFormatter) const {
    std::cout << s << "InertiaRatiosFactor(" << keyFormatter(this->key())
              << "  measured: " << prior_.vector() << std::endl;
    this->noiseModel_->print("  noise model: ");
  }

  /** equals */
  virtual bool equals(const NonlinearFactor& expected,
                      double tol = 1e-9) const {
    const This* e = dynamic_cast<const This*>(&expected);
    return e != NULL && Base::equals(*e, tol) &&
           this->prior_.equals(e->prior_, tol);
    // return e != NULL && Base::equals(*e, tol) &&
    //        traits<Point3>::Equals(this->prior_, e->prior_, tol);
  }

  /** implement functions needed to derive from Factor */

  /** vector of errors */
  Vector evaluateError(const InertiaRatios& ir,
                       boost::optional<Matrix&> H1 = boost::none) const {
    return ir.vector() - prior_.vector();
  }

  /** return the prior */
  const InertiaRatios& prior() const { return prior_; }

 private:
  /** Serialization function */
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
    ar& BOOST_SERIALIZATION_NVP(prior_);
  }
};

}  // namespace gtsam
