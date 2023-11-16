/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/5/6
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_TETMESHFACTORY_H_
#define CVT_INCLUDE_CVT_TETMESHFACTORY_H_

#include "common.h"
#include <vector>
#include <string>

namespace cvt{
namespace tet {
class TetMeshGenerator {
public:
    virtual VF operator()(int cid, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) const = 0;
};

class TetWildMeshGenerator : public TetMeshGenerator {
    double tet_radio;
public:
    explicit TetWildMeshGenerator(double tet_radio);
    VF operator()(int cid, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) const override;
};

class NotGenerateTetMesh: public TetMeshGenerator{
public:
    VF operator()(int cid, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) const override {
        return VF();
    }
};

class TetMeshLoader : public TetMeshGenerator {
    const std::string prefix;
public:
    explicit TetMeshLoader(const std::string &prefix);
    VF operator()(int cid, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) const override;
};

void SaveTetMesh(const std::string &prefix, const std::vector<VF> &tetmesh);

}
}


#endif //CVT_INCLUDE_CVT_TETMESHFACTORY_H_
