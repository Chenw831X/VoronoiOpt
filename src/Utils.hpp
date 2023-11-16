//
// Created by Wei Chen on 3/2/22
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include <Eigen/Eigen>

namespace SIM{

class Utils{
public:
    // get Legendre Polynomial
    static double getLegendrePolynomial(double x, int p);
    // get Legendre Polynomial Derivative
    static double getLegendrePolynomialDerivative(double x, int p);
    // evaluate integrated Legendre basis function
    static double evalPhi(int p, double x);
    // evaluate derivative of integrated Legendre basis function
    static double evalPhiDerivative(int p, double x);

    // get Gauss Quadrature Coordinates
    static void getGaussQuadratureCoordinates(int GO, Eigen::VectorXd &GP);
    // get Gauss Quadrature Weights
    static void getGaussQuadratureWeights(int GO, Eigen::VectorXd &w);

    // evaluate Cubic Bernstein basis, where 0 <= t <= 1
    static void evalCubicBernstein(double t, Eigen::RowVector4d &phi);

    // compute the barycentric weights for a point p in an n-side polygon
    // write answer to weight(rowI, :)
    static void compute_barycentric(const Eigen::MatrixXd &polygon, const Eigen::RowVector3d &p,
                                    int rowI, Eigen::MatrixXd &weight);

    static double cotangent(const Eigen::RowVector3d &a, const Eigen::RowVector3d &b,
                            const Eigen::RowVector3d &c);

    static bool writeMatrixXd(const std::string &filePath, const Eigen::MatrixXd &mat);
    static bool writeMatrixXi(const std::string &filePath, const Eigen::MatrixXi &mat);
    static bool readMatrix(const std::string &filePath, int row, int col, Eigen::MatrixXd &mat);

    // write mesh information to '.txt', prepare for meshio
    // V.size(): number of elements in mesh
    // m = 1: V[i] is (1, 3), point mesh
    // m = 2: V[i] is (2, 3), line mesh
    // m = 4: V[i] is (4, 3), quad mesh
    // m = 8: V[i] is (8, 3), hex mesh
    static bool writeMesh(const std::string &filePath, const std::vector<Eigen::MatrixXd> &V,
                          int m);

    static bool writeMesh0(const std::string &filePath,
                           const std::vector<std::vector<Eigen::MatrixXd>> &polygon,
                           const std::vector<Eigen::MatrixXd> &micV,
                           const std::vector<Eigen::MatrixXi> &micT);

    static bool readMesh(const std::string &filePath,
                         std::vector<std::vector<Eigen::MatrixXd>> &polygon,
                         Eigen::MatrixXd &node,
                         std::vector<Eigen::VectorXi> &eDof,
                         std::vector<std::vector<Eigen::VectorXi>> &nid,
                         std::vector<std::vector<std::vector<Eigen::VectorXi>>> &dup,
                         std::vector<Eigen::MatrixXd> &micV,
                         std::vector<Eigen::MatrixXi> &micT,
                         std::vector<std::vector<Eigen::VectorXi>> &bnid,
                         std::vector<std::vector<Eigen::VectorXi>> &bnid_deDup,
                         std::vector<std::vector<Eigen::VectorXi>> &bnid_deDup_ID,
                         std::vector<int> &bnNode,
                         std::vector<Eigen::VectorXi> &reorderVec,
                         int &nEle, int &nNode, std::vector<int> &nFacePerPoly,
                         std::vector<int> &eleDofNum, std::vector<int> &mic_nV,
                         std::vector<int> &mic_nT);

    static void writeTriVTK(const std::string &path, const Eigen::MatrixXd &V, const Eigen::MatrixXi &T,
                        const std::vector<double> &cell_data = {}, const std::vector<double> &v_data = {});

};

} // namespace SIM

#endif // UTILS_HPP