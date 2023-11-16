//
// Created by Wei Chen on 3/2/22
//

#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <Eigen/Eigen>

namespace SIM{

class Config{
public:
    int example;

    // triangle mesh file path of physical domain
    std::string workspace;
    std::string PD_FilePath;
    std::string cvt_mesh;
    std::string seeds_path;
    std::string macro_prefix;
    int macro_size;
    std::string tet_prefix;
    // polygon mesh file path
    std::string meshFile_;
    std::string meshFile;

    // Material Setup
    double YM1;
    double YM0;
    double PR;

    // penalty parameter in YM
    double penaltyYM;

    // Boundary Condition Setup
    std::vector<Eigen::Matrix<double, 2, 3>> NBCRelBBox;
    std::vector<Eigen::Vector3d> NBCVal;

    std::vector<Eigen::Matrix<double, 2, 3>> DBCRelBBox;
    std::vector<Eigen::Vector3d> DBCVal;

    // penalty parameter of Penalty method in DBC
    double penaltyDBC;

    double t, tmin, tmax;
    double scalarE;
    double volfrac;

public:
    Config();
};

} // namespace SIM

#endif // CONFIG_HPP