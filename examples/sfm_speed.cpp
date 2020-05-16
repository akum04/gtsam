/* ----------------------------------------------------------------------------

 * GTSAM Copyright 2010, Georgia Tech Research Corporation,
 * Atlanta, Georgia 30332-0415
 * All Rights Reserved
 * Authors: Frank Dellaert, et al. (see THANKS for the full author list)

 * See LICENSE for the license information

 * -------------------------------------------------------------------------- */

/**
 * @file    SFMExampleExpressions.cpp
 * @brief   A structure-from-motion example done with Expressions
 * @author  Frank Dellaert
 * @author  Duy-Nguyen Ta
 * @date    October 1, 2014
 */

/**
 * This is the Expression version of SFMExample
 * See detailed description of headers there, this focuses on explaining the AD part
 */

// The two new headers that allow using our Automatic Differentiation Expression framework
#include <gtsam/slam/expressions.h>
#include <gtsam/nonlinear/ExpressionFactorGraph.h>

// Header order is close to far
// #include <examples/SFMdata.h>
#include <gtsam/geometry/Point2.h>
#include <gtsam/nonlinear/DoglegOptimizer.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/slam/dataset.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/ProjectionFactor.h>
#include <gtsam/geometry/SimpleCamera.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;
using namespace gtsam;
using namespace noiseModel;

std::vector<gtsam::Point3> createSatPoints() {

  // Create the set of ground-truth landmarks
  std::vector<gtsam::Point3> points;

  string points_loc = findExampleDataFile("speed_points.txt");
  points.push_back(Point3(0.0, 0.0, 0.0));
  
  ifstream points_file(points_loc.c_str());
  cout << "Reading satellite points" << endl;
  int point_id;
  MatrixRowMajor m(3,1);
  //read camera pose parameters and use to make initial estimates of camera poses
  while (points_file >> point_id) {
    for (int i = 0; i < 3; i++) {
      points_file >> m.data()[i];
    } 
    points.push_back(Point3(m(0), m(1), m(2)));
  }
  return points;
}

/* ************************************************************************* */
int main(int argc, char* argv[]) {

  // Define the camera calibration parameters
  Cal3_S2::shared_ptr K(new Cal3_S2(3003.4, 3003.4, 0.0, 960.0, 600.0));

  // Define the camera observation noise model
  noiseModel::Isotropic::shared_ptr measurementNoise = noiseModel::Isotropic::Sigma(2, 10.0); // one pixel in u and v

  /* ************************************************************************* */



  // Create a factor graph
  NonlinearFactorGraph graph;

  // Add a prior on pose x1. This indirectly specifies where the origin is.
  noiseModel::Diagonal::shared_ptr poseNoise = noiseModel::Diagonal::Sigmas((Vector(6) << Vector3::Constant(0.05), Vector3::Constant(0.01)).finished()); // 5cm std on x,y,z 0.01 rad on roll,pitch,yaw

  std::vector<gtsam::Pose3> poses;

  // Create the set of ground-truth landmarks
  vector<Point3> points = createSatPoints();

  Values initial_estimate;
  string pose_loc = findExampleDataFile("speed_camera_poses.txt");

  ifstream pose_file(pose_loc.c_str());
  cout << "Reading camera poses" << endl;
  int pose_id;
  MatrixRowMajor m(7,1);
  Pose3 tmpPose;
  //read camera pose parameters and use to make initial estimates of camera poses
  while (pose_file >> pose_id) {
    for (int i = 0; i < 7; i++) {
      pose_file >> m.data()[i];
    } 
    tmpPose = Pose3(Rot3::Quaternion(m(0), m(1), m(2), m(3)), Point3(m(4), m(5), m(6)));
    cout << tmpPose.inverse().rotation().quaternion() << tmpPose.inverse().translation() << endl;
    poses.push_back(tmpPose.inverse());
    graph.emplace_shared<PriorFactor<Pose3> >(Symbol('x', pose_id), tmpPose.inverse() , poseNoise); // add directly to graph

    // Pose3 init = Pose3(Rot3::Quaternion(1.0, 0.0, 0.0, 0.0), Point3(0.0, 0.0, 0.0));
    // SimpleCamera camera(tmpPose.inverse(), *K);
    // Point2 pose_measurement = camera.project(init.translation());
    // cout << pose_measurement << endl;
  }


  string lmk_loc = findExampleDataFile("speed_landmark.txt");

  ifstream lmk_file(lmk_loc.c_str());
  cout << "Reading camera poses" << endl;
  int lmk_id;
  MatrixRowMajor lm(4,1);
  //read camera pose parameters and use to make initial estimates of camera poses
  while (lmk_file >> lmk_id) {
    for (int i = 0; i < 4; i++) {
      lmk_file >> lm.data()[i];
    }
    Point2 measurement = Point2(lm(2), lm(3));
    graph.emplace_shared<GenericProjectionFactor<Pose3, Point3, Cal3_S2> >(measurement, measurementNoise, Symbol('x', lm(0)), Symbol('l', lm(1)), K);
  }

  // Because the structure-from-motion problem has a scale ambiguity, the problem is still under-constrained
  // Here we add a prior on the position of the first landmark. This fixes the scale by indicating the distance
  // between the first camera and the first landmark. All other landmark positions are interpreted using this scale.
  noiseModel::Isotropic::shared_ptr pointNoise = noiseModel::Isotropic::Sigma(3, 0.01);
  graph.emplace_shared<PriorFactor<Point3> >(Symbol('l', 0), Point3(0.0, 0.0, 0.0), pointNoise); // add directly to graph
  graph.print("Factor Graph:\n");

  // Create the data structure to hold the initial estimate to the solution
  // Intentionally initialize the variables off from the ground truth
  Values initialEstimate;
  for (size_t i = 0; i < poses.size(); ++i)
    initialEstimate.insert(Symbol('x', i), poses[i].compose(Pose3(Rot3::Rodrigues(-0.01, 0.01, 0.01), Point3(0.05, -0.05, 0.05))));
  for (size_t j = 0; j < points.size(); ++j){
        cout << j << endl;
        initialEstimate.insert<Point3>(Symbol('l', j), points[j] + Point3(-0.1, 0.1, 0.1));
  }
  initialEstimate.print("Initial Estimates:\n");

  /* Optimize the graph and print results */
  Values result = DoglegOptimizer(graph, initialEstimate).optimize();
  result.print("Final results:\n");
  cout << "initial error = " << graph.error(initialEstimate) << endl;
  cout << "final error = " << graph.error(result) << endl;

  return 0;
}
/* ************************************************************************* */