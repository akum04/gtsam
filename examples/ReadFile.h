#include <Eigen/Dense>
#include <fstream>
#include <vector>

#include <gtsam/geometry/Point3.h>
#include <gtsam/geometry/Pose3.h>

using namespace Eigen;

template <typename M>
M load_csv(const std::string& path) {
  std::ifstream indata;
  indata.open(path);
  std::string line;
  std::vector<double> values;
  uint rows = 0;
  while (std::getline(indata, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    while (std::getline(lineStream, cell, ',')) {
      values.push_back(std::stod(cell));
    }
    ++rows;
  }
  return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime,
                          M::ColsAtCompileTime, RowMajor>>(
      values.data(), rows, values.size() / rows);
}

/* ************************************************************************* */
std::vector<gtsam::Point3> createPoints1() {
  // Create the set of ground-truth landmarks
  std::vector<gtsam::Point3> points;
  points.push_back(gtsam::Point3(1.0, 1.0, 1.0));
  points.push_back(gtsam::Point3(-1.0, 1.0, 1.0));
  points.push_back(gtsam::Point3(-1.0, -1.0, 1.0));
  points.push_back(gtsam::Point3(1.0, -1.0, 1.0));
  points.push_back(gtsam::Point3(1.0, 1.0, -1.0));
  points.push_back(gtsam::Point3(-1.0, 1.0, -1.0));
  points.push_back(gtsam::Point3(-1.0, -1.0, -1.0));
  points.push_back(gtsam::Point3(1.0, -1.0, -1.0));

  return points;
}
/* ************************************************************************* */
// std::vector<gtsam::Pose3> createChaserPoses(
//             Eigen::MatrixXd& A,
//             gtsam::Rot3::Ypr& deltaRot,
//             int steps = 8) {

//   // Create the set of ground-truth poses
//   // Default values give a circular trajectory, radius 30 at pi/4 intervals,
//   always facing the circle center Eigen::Vector3d initPoint; initPoint <<
//   A(1,1), A(1,2), A(1,3) ; const gtsam::Pose3& init =
//   gtsam::Pose3(gtsam::Rot3::Ypr(M_PI/2,0,-M_PI/2), gtsam::Point3(initPoint
//   ));

//   std::vector<gtsam::Pose3> poses;
//   int i = 1;
//   Eigen::Vector3d deltaPoint;
//   poses.push_back(init);
//   for(; i < steps; ++i) {
//     deltaPoint << A(i,1), A(i,2), A(i,3) ;
//     const gtsam::Pose3& delta = gtsam::Pose3(deltaRot,
//     gtsam::Point3(deltaPoint)); poses.push_back(poses[i-1].compose(delta));
//   }

//   return poses;
// }
