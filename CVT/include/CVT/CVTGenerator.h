/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/4/17
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT__CVTGENERATOR_H_
#define CVT__CVTGENERATOR_H_
#include <geogram1/geogram/delaunay/periodic_delaunay_3d.h>
#include <geogram1/geogram/delaunay/delaunay_3d.h>
#include "VoronoiCell.h"
#include <Eigen/Geometry>
#include <tuple>
#include "common.h"

struct Vec4Cmp {
    bool operator()(const Eigen::Vector4d &a, const Eigen::Vector4d &b) const {
        if ((a - b).norm() < 1e-6) return false;

        return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y()) ||
            (a.x() == b.x() && a.y() == b.y() && a.z() < b.z()) ||
            (a.x() == b.x() && a.y() == b.y() && a.z() == b.z() && a.w() < b.w());
    };
};

class CVTGenerator {
public:
    GEO::PeriodicDelaunay3d *delaunay_;
    GEO::PeriodicDelaunay3d::IncidentTetrahedra W_;
    std::vector<double> points_;
    std::set<Eigen::Vector4d, Vec4Cmp> convex_planes;

    Eigen::MatrixXd V;
    const Eigen::MatrixXi &F;
    Eigen::AlignedBox<double, 3> bbox;

    Eigen::Vector3d transformFactor = Eigen::Vector3d(0, 0, 0);
    Eigen::Vector3d scaleFactor = Eigen::Vector3d(1, 1, 1);
    double scale = 1;

    explicit CVTGenerator(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
    explicit CVTGenerator(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::AlignedBox3d box);

    std::vector<cvt::VoronoiCell> get_patched_cvt(int npoints, int loyd_iter_n = 100, double split_angle = 0);

    std::vector<cvt::VoronoiCell> get_patched_cvt(const std::vector<Eigen::Vector3d> &points,
                                                  int loyd_iter_n = 100,
                                                  bool cut = true,
                                                  bool split = false,
                                                  double split_angle = 0,
                                                  bool return_all = false,
                                                  bool never_bool = false);

    /**
     * 细节更多的cvt
     * @param points
     * @return
     */
    std::vector<std::pair<cvt::VoronoiCell,
                          cvt::CVTExtraData>> get_detailed_cvt(const std::vector<Eigen::Vector3d> &points);

    void setup_convex_hull();

    void inject_convex_planes(std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> triangles);
    void inject_convex_planes(std::vector<Eigen::Vector4d> planes);

    std::vector<cvt::VoronoiCell> get_cvt(int npoints, int loyd_iter_n = 100);
    std::vector<cvt::VoronoiCell> get_cvt(const std::vector<GEO::vec3> &points, int loyd_iter_n = 100);

    void get_cell(GEO::index_t v, GEO::ConvexCell &C);

    void Lloyd_iteration();

    cvt::VoronoiCell parse_cell(GEO::index_t v);
};

#endif //CVT__CVTGENERATOR_H_
