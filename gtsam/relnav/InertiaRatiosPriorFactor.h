/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  InertiaRatiosPriorFactor.h
 *  @author Arunkumar Rathinam
 **/
#pragma once

#include <gtsam/base/Testable.h>
#include <gtsam/base/VectorSpace.h>

#include <gtsam/nonlinear/NonlinearFactor.h>
#include "InertiaRatios.h"

namespace gtsam {

/**
 * A class for a soft prior on any Value type
 * @addtogroup SLAM
 */
template <class InertiaRatios>
class InertiaRatiosPriorFactor : public NoiseModelFactor1<InertiaRatios> {
 private:
  typedef InertiaRatiosPriorFactor<InertiaRatios> This;
  typedef NoiseModelFactor1<InertiaRatios> Base;

  InertiaRatios prior_; /** The measurement */

  /** concept check by type */
  GTSAM_CONCEPT_TESTABLE_TYPE(InertiaRatios)
  GTSAM_CONCEPT_POSE_TYPE(InertiaRatios)
 public:
  /// shorthand for a smart pointer to a factor
  typedef typename boost::shared_ptr<InertiaRatiosPriorFactor<InertiaRatios> >
      shared_ptr;

  /** default constructor - only use for serialization */
  InertiaRatiosPriorFactor() {}

  virtual ~InertiaRatiosPriorFactor() {}

  /** Constructor */
  InertiaRatiosPriorFactor(Key key, const InertiaRatios& prior,
                           const SharedNoiseModel& model)
      : Base(model, key), prior_(prior) {}

  /// @return a deep copy of this factor
  virtual gtsam::NonlinearFactor::shared_ptr clone() const {
    return boost::static_pointer_cast<gtsam::NonlinearFactor>(
        gtsam::NonlinearFactor::shared_ptr(new This(*this)));
  }

  /** implement functions needed for Testable */

  /** print */
  virtual void print(const std::string& s, const KeyFormatter& keyFormatter =
                                               DefaultKeyFormatter) const {
    std::cout << s << "PriorFactor on " << keyFormatter(this->key()) << "\n";
    prior_.print("  prior mean: ");
    this->noiseModel_->print("  noise model: ");
  }

  /** equals */
  virtual bool equals(const NonlinearFactor& expected,
                      double tol = 1e-9) const {
    const This* e = dynamic_cast<const This*>(&expected);
    return e != NULL && Base::equals(*e, tol) &&
           this->prior_.equals(e->prior_, tol);
  }

  /** implement functions needed to derive from Factor */

  /** vector of errors */
  Vector evaluateError(const InertiaRatios& p,
                       boost::optional<Matrix&> H = boost::none) const {
    return ir.vector() - prior_.vector();
  }

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
