//
// Cube.hpp
// fcm_linear_cpp
//
// Created by Wei Chen on 12/2/21
//

#ifndef CUBE_HPP
#define CUBE_HPP

#include <Eigen/Eigen>

namespace SIM{

// A Cube has a Vertex distribution as follows:
//     3        2
//     o--------o
//    /|       /|
// 7 / |    6 / |
//  o--------o  |
//  |0 o-----|--o 1
//  | /      | /
//  |/       |/
//  o--------o
//  4        5
// 01 is x-direction, 03 is y-direction, 04 is z-direction
class Cube{
public:
    Eigen::MatrixXd Coords;
    
public:
    Eigen::Matrix3d jac; // jacobian matrix
    Eigen::Matrix3d invJac; // inversion of jacobian matrix
    double detJac; // determinant of jacobian matrix
    Eigen::RowVector3d center; // coordinates of element center

public:
    Cube(const Eigen::MatrixXd &p_Coords);

public:
    void debug(void);

    // map local coordinates (standard hex) to global coordinates
    void mapLocalToGlobal(const Eigen::RowVector3d &local, Eigen::RowVector3d &global) const;
    // map global coordinates to local coordiantes (standard hex)
    void mapGlobalToLocal(const Eigen::RowVector3d &global, Eigen::RowVector3d &local) const;
};

} //namespace SIM

#endif // CUBE_HPP