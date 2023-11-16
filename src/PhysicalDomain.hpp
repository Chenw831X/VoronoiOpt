//
// Created by Wei Chen on 3/2/22
//

#ifndef PHYSICALDOMAIN_HPP
#define PHYSICALDOMAIN_HPP

#include <Eigen/Eigen>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/fast_winding_number.h>

#include "Triangle.hpp"

namespace SIM{

class PhysicalDomain{
private:
    std::string filePath;

    igl::FastWindingNumberBVH fwn_bvh;

public:
    // Boundary Conditions
    std::vector<Eigen::Matrix<double, 2, 3>> NBCRelBBox;
    std::vector<Eigen::Matrix<double, 2, 3>> NBCBBox;
    std::vector<Eigen::Vector3d> NBCVal;

    std::vector<Eigen::Matrix<double, 2, 3>> DBCRelBBox;
    std::vector<Eigen::Matrix<double, 2, 3>> DBCBBox;
    std::vector<Eigen::Vector3d> DBCVal;

    int nNBC;
    std::vector<std::pair<Eigen::RowVector3d, int>> NBC;
    int nDBC;
    std::vector<std::pair<int, int>> DBC;

public:
    Eigen::MatrixXd V;
    Eigen::MatrixXd V1;
    Eigen::MatrixXi F;
    std::vector<Triangle> Tri;
    int numV, numF;

    Eigen::Matrix<double, 2, 3> bbox; // bounding box of physical domain
    Eigen::RowVector3d lenBBox; // length of bounding box in 3 dimensions

public:
    PhysicalDomain(const std::string &p_filePath,
        const std::vector<Eigen::Matrix<double, 2, 3>> &p_NBCRelBBox,
        const std::vector<Eigen::Vector3d> &p_NBCVal,
        const std::vector<Eigen::Matrix<double, 2, 3>> &p_DBCRelBBox,
        const std::vector<Eigen::Vector3d> &p_DBCVal);
    
public:
    // after modification of V, update bbox, triangle and fwn_bvh
    void update();

    // given query Q (global coordinates), return W (their domain index)
    // 0: fictitious domain
    // 1: physical domain
    // domain index corresponds to the material index
    void getDomainID(const Eigen::MatrixXd &Q, Eigen::VectorXi &domainID) const;
    // initialize boundary conditions
    void initBC();

    // write triangle mesh based on current V and F
    void writeOBJ(int loop) const;
    // debug
    void writeNBC();
    void writeDBC_obj();
};

} // namespace SIM

#endif // PHYSICALDOMAIN_HPP