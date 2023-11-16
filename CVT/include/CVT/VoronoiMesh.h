/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/4/29
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_VORONOIMESH_H_
#define CVT_INCLUDE_CVT_VORONOIMESH_H_

#include "VoronoiCell.h"
#include <tuple>
#include <Eigen/Geometry>

namespace cvt {


struct VoronoiCellEdgeWrapper {
    static double eps;
    std::vector<Eigen::Vector3d> path;    // 途径的顶点

    PntMat V;
    FaceMat F;

    std::vector<int> cells;               // 被共享的cell

public:
#if 0
    bool equals(const std::vector<Eigen::Vector3d> &other) {
        if (((path.front() - other.front()).norm() < eps) && ((path.back() - other.back()).norm() < eps) ||
            ((path.front() - other.back()).norm() < eps) && ((path.back() - other.front()).norm() < eps)){
            // iff two points coincides
            throw std::runtime_error("Unimpl");
        } else{
            return false;
        }
//        if(path.size() != other.size()) return false;
//
//        bool eq = true;
//
//        for (int i = 0; i < path.size(); ++i) {
//            auto &v1 = path[i];
//            auto &v2 = other[i];
//            if((v1 - v2).norm() >= eps){
//                eq = false;
//                break;
//            }
//        }
//
//        if(eq) return true;
//
//        for (int i = 0; i < path2.size(); ++i) {
//            auto &v1 = path2[i];
//            auto &v2 = other[i];
//            if((v1 - v2).norm() >= eps) return false;
//        }
//
//        return true;
    }
#endif

    bool equals(const std::vector<Eigen::Vector3d> &other, const PntMat &otherV, const FaceMat &otherF);

public:
    Eigen::AlignedBox3d bbox;
};

using cell_edge_ptr = std::shared_ptr<VoronoiCellEdgeWrapper>;

/**
 * Cell connectivity
 */

using RodVec = std::tuple<Eigen::Vector3d, double, std::pair<Eigen::Vector3d, Eigen::Vector3d>, double, Eigen::Matrix3d>;  // <center, L, <a, b>, t, Rotxyz>
struct VoronoiMesh;

struct VoronoiCellWrapper{
    VoronoiCell cell;

    std::vector<cell_edge_ptr> edges;

    std::vector<RodVec> voro_edge_vectors(const VoronoiMesh &mesh) const;
};


struct VoronoiMesh {
    std::vector<VoronoiCellWrapper> cells_connectivities;

    std::vector<cell_edge_ptr> cells_edges;
public:
    explicit VoronoiMesh(const std::vector<VoronoiCell> &cells);

    explicit VoronoiMesh() = default;

private:
    void build_connectivity();

    void build_patch_connectivity();

public:
    std::vector<RodVec> voro_edge_vectors() const;
    std::vector<RodVec> voro_detail_edge_vectors() const;
};
}

#endif //CVT_INCLUDE_CVT_VORONOIMESH_H_
