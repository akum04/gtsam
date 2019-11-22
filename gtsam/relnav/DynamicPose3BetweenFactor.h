
/**
 * @file DynamicPose3BetweenFactor.h
 * @brief
 * @author Arunkumar Rathinam
 */

#pragma once

#include <gtsam/geometry/Cal3_S2.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam_unstable/geometry/InvDepthCamera3.h>

namespace gtsam {

/**
 * Ternary factor representing a visual measurement that includes inverse depth
 */
template <class DynamicPose3, class InertiaRatios>
class DynamicPose3BetweenFactor
    : public NoiseModelFactor3<DynamicPose3, DynamicPose3, InertiaRatios> {
 protected:
  // Keep a copy of measurement and calibration for I/O
  DynamicPose3 measured1_, measured2_;  ///< 2D measurement
  InertiaRatios ir_;
  double dt;

 public:
  /// shorthand for base class type
  typedef NoiseModelFactor3<DynamicPose3, DynamicPose3, InertiaRatios> Base;

  /// shorthand for this class
  typedef DynamicPose3BetweenFactor<DynamicPose3, DynamicPose3, InertiaRatios>
      This;

  /// shorthand for a smart pointer to a factor
  typedef boost::shared_ptr<This> shared_ptr;

  /// Default constructor
  DynamicPose3BetweenFactor()
      : measured_()) {}

  /**
   * Constructor
   * TODO: Mark argument order standard (keys, measurement, parameters)
   * @param measured is the 2 dimensional location of point in image (the
   * measurement)
   * @param model is the standard deviation
   * @param pose1Key is the index of the camera DynamicPose3
   * @param pose2Key is the index of the DynamicPose3
   * @param irKey is the index of inverse depth
   * @param K shared pointer to the constant calibration
   */
  DynamicPose3BetweenFactor(const DynamicPose3& p1measured,
                            const DynamicPose3& p2measured,
                            const InertiaRatios& irmeasured,
                            const SharedNoiseModel& model, const Key pose1Key,
                            Key pose2Key, Key irKey, double timestep)
      : Base(model, pose1Key, pose2Key, irKey),
        measured1_(p1measured),
        measured2_(p2measured),
        ir_(irmeasured),
        dt(timestep) {}

  /** Virtual destructor */
  virtual ~DynamicPose3BetweenFactor() {}

  /**
   * print
   * @param s optional string naming the factor
   * @param keyFormatter optional formatter useful for printing Symbols
   */
  void print(const std::string& s = "DynamicPose3BetweenFactor",
             const KeyFormatter& keyFormatter = DefaultKeyFormatter) const {
    Base::print(s, keyFormatter);
  }

  /// equals
  virtual bool equals(const NonlinearFactor& p, double tol = 1e-9) const {
    const This* e = dynamic_cast<const This*>(&p);
    return e && Base::equals(p, tol);
    // &&
    //        traits<Point2>::Equals(this->measured_, e->measured_, tol) &&
    //        this->K_->equals(*e->K_, tol);
  }

  /// Evaluate error h(x)-z and optionally derivatives
  Vector evaluateError(const DynamicPose3& p1measured,
                       const DynamicPose3& p2measured,
                       const InertiaRatios& irmeasured,
                       boost::optional<Matrix&> H1 = boost::none,
                       boost::optional<Matrix&> H2 = boost::none,
                       boost::optional<Matrix&> H3 = boost::none) const {
    Eigen::VectorXd err =
        p2measured.computeStateChange(p1measured, dt, irmeasured);
    return err;
  }

 private:
  /// Serialization function
  friend class boost::serialization::access;
  template <class ARCHIVE>
  void serialize(ARCHIVE& ar, const unsigned int /*version*/) {
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(Base);
    ar& BOOST_SERIALIZATION_NVP(measured1_);
    ar& BOOST_SERIALIZATION_NVP(measured2_);
    ar& BOOST_SERIALIZATION_NVP(ir_);
    ar& BOOST_SERIALIZATION_NVP(dt);
  }
};
}  // namespace gtsam
