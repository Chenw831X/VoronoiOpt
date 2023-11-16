/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/5/2
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_CVTGEOMETRY_H_
#define CVT_INCLUDE_CVT_CVTGEOMETRY_H_

#include <tuple>
#include "VoronoiCell.h"

namespace cvt {
// Voronoi, Tet V, Tet Mesh, Tet cell value, surface v
using BackgroundCellType = std::tuple<cvt::VoronoiCell,   // 0: Voronoi
                                      Eigen::MatrixXd,    // 1: Tet V
                                      Eigen::MatrixXi,    // 2: Tet Mesh
                                      Eigen::MatrixXd,    // 3: Tet cell value
                                      Eigen::MatrixXi>;   // 4: surface v

class CVTGeometry {
    const Eigen::MatrixXd &V;
    const Eigen::MatrixXi &F;

    std::vector<BackgroundCellType> background_cells;

public:
    CVTGeometry(const Eigen::MatrixXd &V,
                const Eigen::MatrixXi &F);

    /**
     *
     * @param V 网格顶点
     * @param F 网格面
     * @param expand_v bounding box 扩大值，
     * @param sample_cnt voronoi采样数量，输出的voronoi种子点数量可能比这个小
     * @param radius     默认杆半径
     * @param llyod_iter_n   cvt优化迭代次数
     * @param tet_radio    生成四面体的分辨率，为一个voronoi单元bounding box对角线的百分比(0~1)
     * @return
     */
    std::vector<BackgroundCellType> init(double expand_v,
                                         int sample_cnt,
                                         double radius,
                                         int llyod_iter_n,
                                         double tet_radio);

    std::vector<Eigen::VectorXd> calc_rho(const Eigen::MatrixXd &X);
};
}

#endif //CVT_INCLUDE_CVT_CVTGEOMETRY_H_
