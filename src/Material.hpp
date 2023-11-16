//
// Created by Wei Chen on 3/2/22
//

#ifndef MATERIAL_HPP
#define MATERIAL_HPP

#include <Eigen/Eigen>

namespace SIM{

class Material{
public:
    double YM, PR; // Young's modulus, Poisson's ratio
    Eigen::Matrix<double, 6, 6> D; // constitutive matrix

public:
    Material(double p_YM, double p_PR);

    // compute element stiffness matrix of 3D linear elastic material
    // using 2-nd Gauss integration
    // @input
    //  a, b, c: half size of element
    //  D      : constitutive matrix, (6, 6)
    // @output
    //  Ke     : element matrix, (24, 24)
    static void computeKe(double a, double b, double c, const Eigen::Matrix<double, 6, 6> &D,
                          Eigen::Matrix<double, 24, 24> &Ke);

    // compute element strain-displacement matrix B in 8 Gauss Points
    // @input
    //  a, b, c: half size of element
    // @output
    //  Be     : element strain-displacement matrix, (6, 24)
    static void computeBe(double a, double b, double c, std::vector<Eigen::Matrix<double, 6, 24>> &Be);

    // compute element shape function matrix in a given local point P
    // @input
    //  P: given local point, in [-1, 1] * [-1, 1] * [-1, 1]
    // @output
    //  N: shape function matrix N, (3, 24)
    static void computeN(const Eigen::RowVector3d &P, Eigen::Matrix<double, 3, 24> &N);

    // compute element strain-displacement matrix B in a given local point P
    // @input
    //  a, b, c: half size of element
    //  P: given local point, in [-1, 1] * [-1, 1] * [-1, 1]
    // @output
    //  B: element strain-displacement matrix, (6, 24)
    static void computeB(double a, double b, double c, const Eigen::RowVector3d &P, Eigen::Matrix<double, 6, 24> &B);

    static void computeN_tet(const Eigen::RowVector3d &P, const Eigen::Matrix<double, 4, 3> &X,
                             Eigen::Matrix<double, 3, 12> &N);
    static void computeB_tet(const Eigen::Matrix<double, 4, 3> &X, Eigen::Matrix<double, 6, 12> &B);
    static void computeKe_tet(const Eigen::Matrix<double, 4, 3> &X, const Eigen::Matrix<double, 6, 6> &D,
                              Eigen::Matrix<double, 12, 12> &Ke, double &Vol);
};

} // namespace SIM

#endif // MATERIAL_HPP