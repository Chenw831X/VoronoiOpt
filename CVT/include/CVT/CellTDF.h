/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/4/29
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_CELLTDF_H_
#define CVT_INCLUDE_CVT_CELLTDF_H_

#include "VoronoiMesh.h"

namespace cvt {
class CellTDF {
public:
    const std::vector<RodVec> rod_vecs;
public:
    CellTDF(const VoronoiMesh& mesh);

    double phi(const Eigen::Vector3d &p) const;

    double dis(const Eigen::Vector3d &p) const;

    static double H(double v, double eps);
};
}

#endif //CVT_INCLUDE_CVT_CELLTDF_H_
