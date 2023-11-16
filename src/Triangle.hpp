//
// Triangle.hpp
// fcm_linear_cpp
//
// Created by Wei Chen on 12/16/21
//

#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include <Eigen/Eigen>

namespace SIM{

class Triangle{
public:
    // use standard rectangle to fit triangle
    Eigen::Matrix<double, 4, 3> Coords;

public:
    explicit Triangle(const Eigen::Matrix3d &p_Coords);

public:
    // update Coords
    void update(const Eigen::Matrix3d &p_Coords);

    // get two base vector
    Eigen::RowVector3d getBase1(const Eigen::RowVector3d &local) const;
    Eigen::RowVector3d getBase2(const Eigen::RowVector3d &local) const;
    // compute determinant of Jacobian matrix
    double getDetJac(const Eigen::RowVector3d &local) const;
    // get normal vector
    Eigen::RowVector3d getNormal(const Eigen::RowVector3d &local) const;
    // map local coordinates (standard rectangle) to global coordinates
    void mapLocalToGlobal(const Eigen::RowVector3d &local, Eigen::RowVector3d &global) const;
};

} // namespace SIM

#endif // TRIANGLE_HPP