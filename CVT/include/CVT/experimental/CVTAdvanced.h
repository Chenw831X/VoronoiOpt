/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/5/11
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_EXPERIMENTAL_CVTADVANCED_H_
#define CVT_INCLUDE_CVT_EXPERIMENTAL_CVTADVANCED_H_

#include "CVT/common.h"
#include "CVT/CVT.h"
#include "CVT/VoronoiMesh.h"
#include "CVT/CVTSeedFactory.h"
#include "CVT/TetMeshFactory.h"

namespace cvt {
class CVTAlgAdvanced : public CVTAlg {
public:
    explicit CVTAlgAdvanced(const Eigen::MatrixXd &V,
                            const Eigen::MatrixXi &F);

    void step_fd(const Eigen::MatrixXd &X) override;

    std::vector<Eigen::VectorXd> calc_rho_fd() override;

    std::vector<Eigen::VectorXd> calc_rho_fd(Eigen::VectorXi &flag,
                                             const int seed_id,
                                             const Eigen::Vector4d &seed_pos,
                                             const Eigen::Vector3d &search_range);

    Eigen::VectorXd calc_rho_fd_in_single_cell(int cell_id, int updated_seed_id = -1);

public:
    Eigen::MatrixXd current_seeds;
    Eigen::VectorXd current_t;
    std::vector<Eigen::VectorXd> current_rhos;

    // speed up v
    std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> seeds_5_distance_temp;
};
}

#endif //CVT_INCLUDE_CVT_EXPERIMENTAL_CVTADVANCED_H_
