/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/5/4
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_SAMPLING_H_
#define CVT_INCLUDE_CVT_SAMPLING_H_

#include <Eigen/Core>
#include <vector>

namespace cvt {
namespace sampling {
std::vector<std::pair<int, Eigen::Vector2d>> sample_on_mesh(const Eigen::MatrixXd &V,
                                                            const Eigen::MatrixXi &F,
                                                            const int n);

const Eigen::MatrixXd get_points_by_samples(const std::vector<std::pair<int, Eigen::Vector2d>> &samples,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F);


const Eigen::MatrixXd sample_in_mesh(const Eigen::MatrixXd &V,
                                     const Eigen::MatrixXi &F,
                                     int n, double iso=0);


std::vector<std::pair<int, Eigen::Vector2d>> sample_on_corners(const Eigen::MatrixXd &V,
                                                               const Eigen::MatrixXi &F,
                                                               const int n);
};
}

#endif //CVT_INCLUDE_CVT_SAMPLING_H_
