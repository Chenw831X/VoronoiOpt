/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/5/4
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_CVT_H_
#define CVT_INCLUDE_CVT_CVT_H_

#include "common.h"
#include "VoronoiMesh.h"
#include "CVTSeedFactory.h"
#include "TetMeshFactory.h"

namespace cvt{

class CVTAlg {
public:
    explicit CVTAlg(const Eigen::MatrixXd &V,
                    const Eigen::MatrixXi &F);

    virtual void background_mesh(
        const seed::SeedGenerator &seed_generator,
        const tet::TetMeshGenerator &tet_mesh_generator,
        double radius);

    virtual void background_mesh(const std::string &macro_grid_prefix,
                                 const int macro_grid_size,
                                 const tet::TetMeshGenerator &tet_mesh_generator,
                                 double radius);

    void set_shell(const bool use_shell, const double shell_radius);

    void extract_background_macro_unit_meta_data();

    void step(const Eigen::MatrixXd &X);

    virtual void step_fd(const Eigen::MatrixXd &X);

    Eigen::MatrixXd get_vars() const;

    void get_geo_centers(Eigen::MatrixXd & centers, Eigen::VectorXi & flag) const;

    std::vector<Eigen::VectorXd> calc_rho();

    Eigen::VectorXd calc_rho(const Eigen::MatrixXd &P);

    virtual std::vector<Eigen::VectorXd> calc_rho_fd();

    std::pair<int, int> point_loc(const Eigen::Vector3d &p) const;

    std::pair<int, int> point_loc_alternative(const Eigen::Vector3d &p, double tol);
    // alternative to point_loc
    std::vector<std::vector<Eigen::Matrix3d>> tets_inv_bases;

public:
    const Eigen::MatrixXd V;
    const Eigen::MatrixXi F;
    Eigen::AlignedBox3d box;
    Eigen::MatrixXd convexhullV;
    Eigen::MatrixXi convexhullF;
    std::vector<std::vector<Eigen::VectorXi>> tet_boundary_nids;

    std::vector<VF> tet_mesh;
    std::vector<Eigen::AlignedBox3d> macro_mesh_aabb;
    std::vector<Eigen::VectorXd> tet_v_sdf;
    cvt::VoronoiMesh voronoi_mesh;
    cvt::VoronoiMesh bg_voronoi_mesh;
//    std::vector<TriMesh> bg_meshes;

    std::vector<Eigen::Vector3d> init_seeds;
    double eps;

    enum class BackgroundType {
        FromVoronoi,
        FromMesh
    };
    BackgroundType bg_type;

    bool use_shell=false;
    double shell_radius;
};
}

#endif //CVT_INCLUDE_CVT_CVT_H_
