/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 *  @file  KinematicPose3PriorFactor.h
 *  @author Arunkumar Rathinam
 **/
#pragma once

#include <gtsam/base/Testable.h>
#include <gtsam/base/VectorSpace.h>

#include <gtsam/nonlinear/NonlinearFactor.h>
#include "KinematicPose3.h"

namespace gtsam {

/**
 * A class for a soft prior on any Value type
 * @addtogroup SLAM
 */
template <class KinematicPose3>
class KinematicPose3PriorFactor : public NoiseModelFactor1<KinematicPose3> {
 private:
  typedef KinematicPose3PriorFactor<KinematicPose3> This;
  typedef NoiseModelFactor1<KinematicPose3> Base;

  KinematicPose3 prior_; /** The measurement */

  /** concept check by type */
  GTSAM_CONCEPT_TESTABLE_TYPE(KinematicPose3)
  GTSAM_CONCEPT_POSE_TYPE(KinematicPose3)
 public:
  /// shorthand for a smart pointer to a factor
  typedef typename boost::shared_ptr<KinematicPose3PriorFactor<KinematicPose3> >
      shared_ptr;

  /** default constructor - only use for serialization */
  KinematicPose3PriorFactor() {}

  virtual ~KinematicPose3PriorFactor() {}

  /** Constructor */
  KinematicPose3PriorFactor(Key key, const KinematicPose3& prior,
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
  Vector evaluateError(const KinematicPose3& p,
                       boost::optional<Matrix&> H = boost::none) const {
    Vector err = p.vector() - prior_.vector();
    Vector4 q1_tot = p.qTotal();
    Vector4 qm_tot = prior_.qTotal();
    Vector4 dq = p.quaternionDivision(q1_tot, qm_tot);
    Vector3 da = p.quaternion2mrp(dq);

    err.segment<3>(3) = da;

    return err;
  }

  const KinematicPose3& prior() const { return prior_; }

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
