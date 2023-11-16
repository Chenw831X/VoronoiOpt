/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/4/17
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT__COMMON_H_
#define CVT__COMMON_H_

#include <Eigen/Core>
#include <OpenMesh/Core/Mesh/DefaultTriMesh.hh>
#include <OpenMesh/Core/Mesh/DefaultPolyMesh.hh>

namespace cvt {
using Pnt3d = Eigen::Vector3d;
using Vec3d = Eigen::Vector3d;
using Vec4d = Eigen::Vector4d;
using PntMat = Eigen::MatrixXd;
using FaceMat = Eigen::MatrixXi;

// openmesh
using TriMesh = OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraitsDouble>;
using PolyMesh = OpenMesh::PolyMesh_ArrayKernelT<OpenMesh::DefaultTraitsDouble>;

struct VF {
    Eigen::MatrixXd TV;
    Eigen::MatrixXi TT;
};

struct CVTExtraData{
    Eigen::MatrixXd cutSurfV;   // 只是表面的cut
    Eigen::MatrixXi cutSurfF;
    bool intersect;         // 体上是否有交，和上面的不一样
    Eigen::MatrixXd intersectV;
    Eigen::MatrixXi intersectF;
};
}

#endif //CVT__COMMON_H_
