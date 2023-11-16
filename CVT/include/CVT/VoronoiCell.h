/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/4/20
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT__VORONOICELL_H_
#define CVT__VORONOICELL_H_

#include <Eigen/Core>
#include "common.h"
#include <vector>
#include <memory>

namespace cvt {
struct VoronoiCell;

struct Patch;

struct Edge;

struct Edge{
    std::vector<int> vs;  // vertices
    int p0=-1, p1=-1;     // patch id
};

using edge_ptr = std::shared_ptr<Edge>;

struct Patch {
    std::vector<int> face_idx;  // face id

    Vec4d pseudo_plane;         // project plane

    std::vector<edge_ptr> edges;  // boundary edges

    std::vector<std::pair<int, int>> edges_in_order;

    bool surface = false;

public:
    std::vector<OpenMesh::SmartVertexHandle> vertices(const TriMesh &mesh) const;

    cvt::Pnt3d pseudo_pos(const TriMesh &mesh, int v) const;

    cvt::Pnt3d pseudo_pos(const cvt::Pnt3d &v) const;   // project v onto pseudo_plane

    enum class PseudoPlaneType {
        Random,
        First3,
        FitPlane
    };
    void init_pseudo_plane(const TriMesh &mesh, const PseudoPlaneType &type);
};


struct VoronoiCell {
    TriMesh mesh;
    std::vector<Patch> patches;
    std::vector<int> face_patches;
    Vec3d center;   // cell center
    std::vector<edge_ptr> cell_edges;
    std::set<int> corners;  // corner vertices
    bool valid = true;

    double t = 1;  //

public:
    void init_patches();

    static std::vector<VoronoiCell> SplitCellWithComponent(const VoronoiCell &cell);

    /**
     * Compact patches in 0..n
     */
    void compact_patch_idx();

    void divide_patches();

    void split_patches_if_unconnected();

    void split_patches_if_acute_angle(double angle);

public:
    void to_obj(const std::string &filename);

public:
    void build_patch_brep();

public:
    void to_tet();
};
}

#endif //CVT__VORONOICELL_H_
