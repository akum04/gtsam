/* This file will simulate the motion of landmarks based on the corners of
 * the rotating cube*/

// For loading the data
#include "ReadFile.h"
#include "SFMdata1.h"
// #include "matplotlibcpp.h"

#include <gtsam/geometry/Point2.h>
#include <gtsam/geometry/SimpleCamera.h>

#include <gtsam/inference/Symbol.h>
// Isam2 to solve the problem incrementally
#include <gtsam/nonlinear/ISAM2.h>

#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>

#include <gtsam/navigation/Scenario.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/ProjectionFactor.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include <fstream>
#include <iostream>

using namespace std;
using namespace gtsam;

using symbol_shorthand::X;

int main(int argc, char *argv[]) {
  // Define Camera calibration parameters
  Cal3_S2::shared_ptr K(new Cal3_S2(50.0, 50.0, 0.0, 50.0, 50.0));

  // Define the camera observation noise model, 1 pixel std
  auto measurementNoise = noiseModel::Isotropic::Sigma(2, 1.0);

  // Create the set of ground truth landmarks
  std::vector<gtsam::Point3> points = createPoints1();

  cout << "*** Points ***" << endl;
  for (size_t i = 0; i < points.size(); ++i) cout << points[i] << endl;

  // int steps = 8;

  // Create the set of ground-truth poses
  std::vector<gtsam::Pose3> poses;

  // Create an iSAM2 object.
  ISAM2Params parameters;
  parameters.relinearizeThreshold = 0.01;
  parameters.relinearizeSkip = 1;
  ISAM2 isam(parameters);

  // Create a Factor Graph and Values to hold the new data
  NonlinearFactorGraph graph;
  Values initialEstimate;

  Pose3 pose_0 = Pose3(Rot3::RzRyRx(0, 0, 0), Point3(0, 0, 0));
  double angular_velocity = M_PI / 20;
  Vector3 angular_velocity_vector(0, 0, -angular_velocity);
  Vector3 linear_velocity_vector(0, 0, 0);
  auto scenario = ConstantTwistScenario(angular_velocity_vector,
                                        linear_velocity_vector, pose_0);
  // Pose3 delta(Rot3::Rodrigues(0.01, 0.01, -0.01), Point3(0.005, -0.005,
  // 0.001));

  std::vector<gtsam::Point3> new_points;
  for (size_t m = 0; m < points.size(); ++m)
    new_points.push_back(Point3(0, 0, 0));

  // Start with a camera on looking at origin
  const Point3 up(0, 0, 1), target(0, 0, 0);
  const Point3 position(5, 5, 5);
  const SimpleCamera camera = SimpleCamera::Lookat(position, target, up);
  const auto init = camera.pose();
  poses.push_back(init);

  std::cout << "cam_pose: " << init << std::endl;

  // static Pose3 left_cam_tran(Rot3::RzRyRx(0.0, 0.0, 0.0),
  //                            Point3(0.025, 0.0, 0.0));
  //
  // static Pose3 right_cam_tran(Rot3::RzRyRx(0.0, 0.0, 0.0),
  //                             Point3(-0.025, 0.0, 0.0));
  //
  // Pose3 left_cam_pose(init * left_cam_tran);
  // Pose3 right_cam_pose(init * right_cam_tran);
  //
  // std::cout << "left_cam_pose: " << left_cam_pose << std::endl;
  // std::cout << "right_cam_pose: " << right_cam_pose << std::endl;

  static auto kPosePrior = noiseModel::Diagonal::Sigmas(
      (Vector(6) << 0.3, 0.3, 0.3, 0.01, 0.01, 0.01).finished());

  static Pose3 kDeltaPose(Rot3::RzRyRx(-0.001, 0.002, 0.005),
                          Point3(0.15, -0.20, 0.05));

  // include pose to the graph
  graph.emplace_shared<PriorFactor<Pose3>>(Symbol('x', 0), poses[0],
                                           kPosePrior);

  // add odometry noise for the measurement between the poses
  static auto rotationNoise = noiseModel::Diagonal::Sigmas(
      (Vector(6) << 0.0001, 0.0001, 0.0001, 0.01, 0.01, 0.01).finished());

  // trying to create a scenario based prediction
  double delta_t = 0.1;
  size_t frames = 3;

  for (size_t i = 0; i < frames; ++i) {
    double t = i * delta_t;
    auto pose_target = scenario.pose(t);

    initialEstimate.insert(Symbol('t', i), pose_target * kDeltaPose);

    initialEstimate.insert(Symbol('x', i), poses[0] * kDeltaPose);
    poses.push_back(pose_target);

    if (i >= 1) {
      graph.emplace_shared<BetweenFactor<Pose3>>(
          Symbol('t', i - 1), Symbol('t', i),
          Pose3(Rot3::RzRyRx(0.005, 0.005, -angular_velocity),
                Point3(0.005, -0.005, 0.001)),
          rotationNoise);
      graph.emplace_shared<BetweenFactor<Pose3>>(
          Symbol('x', i - 1), Symbol('x', i),
          Pose3(Rot3::RzRyRx(0.001, 0.001, 0.001),
                Point3(0.001, -0.001, 0.001)),
          rotationNoise);
    }

    for (size_t m = 0; m < points.size(); ++m) {
      if (i == 0) {
        initialEstimate.insert(Symbol('l', m), points[m]);
      }
      new_points[m] = pose_target.transformTo(points[m]);
      SimpleCamera camera(poses[0], *K);
      Point2 measurement = camera.project(new_points[m]);
      graph.emplace_shared<GenericProjectionFactor<Pose3, Point3, Cal3_S2>>(
          measurement, measurementNoise, Symbol('x', i), Symbol('l', m), K);
    }
  }
  // initialEstimate.print("Initial Estimates:\n");
  // graph.print("Factor Graph:\n");

  std::cout << "initial error: " << graph.error(initialEstimate) << std::endl;
  Values result =
      LevenbergMarquardtOptimizer(graph, initialEstimate).optimize();
  result.print("Final Result:\n");
  std::cout << "final error: " << graph.error(result) << std::endl;
}
