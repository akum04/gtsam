// \callgraph
#pragma once

#include <gtsam/config.h>

#include <gtsam/base/Lie.h>
#include <gtsam/geometry/BearingRange.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/relnav/RigidBodyDynamics.h>

namespace gtsam {

class GTSAM_EXPORT DynamicPose : public LieGroup<DynamicPose, 12> {
 private:
  RigidBodyDynamics rbd;

 public:
};

}  // namespace gtsam