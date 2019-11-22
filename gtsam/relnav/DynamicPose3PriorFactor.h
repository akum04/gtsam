/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  DynamicPose3PriorFactor.h
 *  @author Arunkumar Rathinam
 **/
#pragma once

#include <gtsam/base/Testable.h>
#include <gtsam/geometry/concepts.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include "DynamicPose3.h"

namespace gtsam {

/**
 * A class for a soft prior on any Value type
 * @addtogroup SLAM
 */
template <class DynamicPose3>
class DynamicPose3PriorFactor : public NoiseModelFactor1<DynamicPose3> {
 private:
  typedef DynamicPose3PriorFactor<DynamicPose3> This;
  typedef NoiseModelFactor1<DynamicPose3> Base;

  DynamicPose3 prior_; /** The measurement */

  /** concept check by type */
  GTSAM_CONCEPT_TESTABLE_TYPE(DynamicPose3)
  GTSAM_CONCEPT_POSE_TYPE(DynamicPose3)
 public:
  /// shorthand for a smart pointer to a factor
  typedef typename boost::shared_ptr<DynamicPose3PriorFactor<DynamicPose3> >
      shared_ptr;

  /** default constructor - only use for serialization */
  DynamicPose3PriorFactor() {}

  virtual ~DynamicPose3PriorFactor() {}

  /** Constructor */
  DynamicPose3PriorFactor(Key key, const DynamicPose3& prior,
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
  Vector evaluateError(const DynamicPose3& p,
                       boost::optional<Matrix&> H = boost::none) const {
    Vector err = p.vectorFull() - prior_.vector();

    Vector4 p1qTot = p.qTotal();
    Vector4 mqTot = prior_.qTotal();
    Vector3d da = p.getMRPdifference(p1qTot, mqTot);
    err.segment<3>(6) = da;

    return err;
    // if (H) (*H) = I_6x6;
    // manifold equivalent of h(x)-z -> log(z,h(x))
    // return prior_.localCoordinates(p);
  }

  const DynamicPose3& prior() const { return prior_; }

 private:
  /** Serialization function */
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
    ar& BOOST_SERIALIZATION_NVP(prior_);
    ar& BOOST_SERIALIZATION_NVP(body_P_sensor_);
  }
};

}  // namespace gtsam
