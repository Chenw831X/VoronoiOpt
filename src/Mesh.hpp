//
// Created by Wei Chen on 4/29/22
//

#ifndef MESH_HPP
#define MESH_HPP

#include <set>
#include <Eigen/Eigen>

#include "PhysicalDomain.hpp"
#include "Material.hpp"
#include "CVT/CVT.h"
#include "CVT/utils.h"
#include "CVT/experimental/CVTAdvanced.h"

#include "EigenLibSolver.hpp"

#ifdef LINSYSSOLVER_USE_CHOLMOD
#include "CHOLMODSolver.hpp"
#endif

namespace SIM {

class Mesh {
public:
    double YM1, YM0;
    double penaltyYM; // penalty parameter of Young's modulus
    double penaltyDBC; // penalty parameter of Penalty method in DBC
    std::shared_ptr<PhysicalDomain> physicalDomain;
    std::shared_ptr<cvt::CVTAlgAdvanced> cvt_alg;
    // number of voronoi cell
    // cell number may not equal element number
    int nCell;
    // number of design variables
    int nVar;

public:
    // segments (n, 6) of each face in each macro element (polygon)
    // n segments are ordered counter clock-wise
    std::vector<std::vector<Eigen::MatrixXd>> polygon;
    // all nodes in polygon mesh, (nNode, 3)
    Eigen::MatrixXd node;
    // macro dofs id in each element
    std::vector<Eigen::VectorXi> eDof;
    // nodes index of each face w.r.t each polygon
    std::vector<std::vector<Eigen::VectorXi>> nid;
    // duplicated index in S-Patches of each node in each face in each polygon
    std::vector<std::vector<std::vector<Eigen::VectorXi>>> dup;

    // V in micro tet mesh in each macro polygon
    std::vector<Eigen::MatrixXd> micV;
    // T in micro tet mesh in each macro polygon
    std::vector<Eigen::MatrixXi> micT;
    // micro nodes index of each face w.r.t tet in each polygon
    std::vector<std::vector<Eigen::VectorXi>> bnid;
    // de-duplicate micro nodes index of each face w.r.t tet in each polygon
    std::vector<std::vector<Eigen::VectorXi>> bnid_deDup;
    // bnid_deDup's index in bnid
    std::vector<std::vector<Eigen::VectorXi>> bnid_deDup_ID;
    // number of micro boundary nodes in each macro polygon
    std::vector<int> bnNode;
    // reorder vector (bnDof, inDof) of micro tet mesh in each macro polygon
    std::vector<Eigen::VectorXi> reorderVec;

public:
    // material
    Material *material = nullptr;

    int nEle; // number of macro elements (polygon)
    int nNode; // number of vertex in this mesh
    int nDof; // number of dofs in this mesh
    // number of face in each macro polygon
    std::vector<int> nFacePerPoly;
    // number of macro dofs in each macro polygon
    std::vector<int> eleDofNum;

    // number of micro nodes in each macro polygon
    std::vector<int> mic_nV;
    // number of micro tets in each macro polygon
    std::vector<int> mic_nT;

    // nodes id in each element
    std::vector<Eigen::VectorXi> eNode;
    // store neighbor nodes id of each node (neighbor: belong to one element)
    std::vector<std::set<int> > vNeighbor;

    // step of finite differential
    double step;

public:
    // transformation matrix for each macro polygon element
    // transformation from micro surface dofs to macro dofs
    std::vector<Eigen::MatrixXd> Phi;

    // element ke of each micro tet in each macro polygon element
    std::vector<std::vector<Eigen::Matrix<double, 12, 12>>> mic_Ke;
    // volume of each micro tet in each macro polygon element
    std::vector<Eigen::VectorXd> mic_Vol;
    // model's initial volume
    double Vol0;

    // dofs id of all micro tets in each macro polygon element
    std::vector<Eigen::VectorXi> mic_eDof;
    // number of micro boundary dofs in each macro polygon element
    std::vector<int> nBDof;
    // number of micro internal dofs in each macro polygon element
    std::vector<int> nIDof;

    std::vector<Eigen::MatrixXd> elementM2;
    // macro element stiffness matrix
    std::vector<Eigen::MatrixXd> elementKe;
    // macro linear system solver
    LinSysSolver<Eigen::VectorXi, Eigen::VectorXd> *linSysSolver = nullptr;

    // macro element load vector
    std::vector<Eigen::VectorXd> elementLoad;
    // load vector
    Eigen::VectorXd load;

    // macro element solution(U) vector
    std::vector<Eigen::VectorXd> elementU;
    // solution(U) vector
    Eigen::VectorXd U;

    // integration order of Boundary Condition (BC) integration
    int integrationOrder_BC = 2;
    Eigen::VectorXd GP_BC; // 1D Gauss Quadrature Points for BC
    Eigen::VectorXd GW_BC; // 1D Gauss Quadrature Weights for BC

#ifdef FCM_NBC
    // micro element index of Gauss Points in NBC
    std::vector<int> NBC_micI;
    // micro N of Gauss Points in NBC
    std::vector<Eigen::Matrix<double, 3, 12>> NBC_micN;
    // integration weight of Gauss Points in NBC
    // std::vector<double> NBC_w;
    // NBC value of Gauss Points in NBC
    std::vector<Eigen::Vector3d> NBC_val;
    // NBC's Gauss points index of each macro element
    std::vector<std::vector<int>> elementNBC;
    bool NBC_flag = false;
#endif // FCM_NBC

#ifdef FCM_DBC
    // micro element index of Gauss Points in DBC
    std::vector<int> DBC_micI;
    // micro N of Gauss Points in DBC
    std::vector<Eigen::Matrix<double, 3, 12>> DBC_micN;
    // micro B of Gauss Points in DBC
    std::vector<Eigen::Matrix<double, 6, 12>> DBC_micB;
    // normal vector Voigt of Gauss Points in DBC
    std::vector<Eigen::Matrix<double, 6, 3>> DBC_DT_mul_normal;
    // integration weight of Gauss Points in DBC
    std::vector<double> DBC_w;
    // DBC value of Gauss Points in DBC
    std::vector<Eigen::Vector3d> DBC_val;
    // DBC's Gauss points index of each macro element
    std::vector<std::vector<int>> elementDBC;
#else // FCM_DBC
    std::vector<int> fixeddofs;
#endif // FCM_DBC

    // number of vertex in physical domain
    int numV;
    // macro element index of vertex in physical domain
    std::vector<int> V_macI;
    // micro element index of vertex in physical domain
    std::vector<int> V_micI;
    // micro N of vertex in physical domain
    std::vector<Eigen::Matrix<double, 3, 12>> V_micN;
    // micro B of vertex in physical domain
    std::vector<Eigen::Matrix<double, 6, 12>> V_micB;
    // macro N of vertex in physical domain
    std::vector<Eigen::MatrixXd> V_macN;
    // macro B of vertex in physical domain
    std::vector<Eigen::MatrixXd> V_macB;
    // vertex index of each macro element
    std::vector<std::vector<int>> elementV;

    // number of vertex in fine tetrahedral mesh
    int fine_num;
    // macro element index of vertex in fine tetrahedral mesh
    std::vector<int> fine_macI;
    // micro element index of vertex in fine tetrahedral mesh
    std::vector<int> fine_micI;
    // micro N of vertex in fine tetrahedral mesh
    std::vector<Eigen::Matrix<double, 3, 12>> fine_micN;
    // macro N of vertex in fine tetrahedral mesh
    std::vector<Eigen::MatrixXd> fine_macN;
    // vertex index of each macro element
    std::vector<std::vector<int>> elementFine;

    // initial fine vertex
    Eigen::MatrixXd fine_V;
    // deformed fine vertex
    Eigen::MatrixXd fine_V_deformed;

    // variables used for locally sensitivity computation
    Eigen::Vector3d fd_bbox;
    std::vector<Eigen::VectorXi> effect_flag;
    std::vector<Eigen::VectorXi> effect_id;
    std::vector<int> effect_num;

public:
    Mesh(double p_YM1, double p_YM0, double p_PR, double p_penaltyYM, double p_penaltyDBC,
         int p_nCell, int p_nVar, std::shared_ptr<PhysicalDomain> p_physicalDomain,
         std::shared_ptr<cvt::CVTAlgAdvanced> p_cvt_alg,
         std::vector<std::vector<Eigen::MatrixXd>> p_polygon,
         std::vector<Eigen::MatrixXd> p_micV,
         std::vector<Eigen::MatrixXi> p_micT
         );

    ~Mesh();

public:
    /**
     * @brief based on polygon, micV, micT, preprocess mesh, i.e.
     * calculate node, eDof, nid, dup, bnid, bnid_deDup, bnid_deDup_ID,
     * bnNode, reorderVec, nEle, nNode, nFacePerPoly, eleDofNum, mic_nV, mic_nT
     */
    void PreprocessMesh();

    void computeFeatures();

    // compute Phi on one face
    void Phi_face(const Eigen::MatrixXd &macV_f, const Eigen::MatrixXd &bV_f, Eigen::MatrixXd &Phi_f);

    // compute Phi
    void computePhi();

    // preprocess NBC, DBC, vertex in physical domain
    void preprocess();

    // prepare for micro integration: compute gp_ID, localGP_K, localGP_B, localGP_w
    void prepareMicroIntegration();

#ifdef USE_CUDA
    // prepare data structure needed by gpu
void prepareGPU();
#endif

    // compute M: transformation from micro surface dofs to micro internal dofs
    // compute and assemble global stiffness matrix, global load vector (Neumann boundary conditions)
    // apply Dirichlet boundary conditions
    void computeSystem(const std::vector<Eigen::VectorXd> &rho);

    void generate_rho(std::vector<Eigen::VectorXd> &rho);
    void solve(); // solve the linear system

    // compute model's current volume: Vol
    double compute_Vol(const std::vector<Eigen::VectorXd> &rho);

    // compute cvt energy
    void cvt_energy(const Eigen::MatrixXd &X, double &E, Eigen::VectorXd &dE);

    // compute dH
    void compute_dH(const Eigen::MatrixXd &X, std::vector<std::vector<Eigen::VectorXd>> &dH);

    // compute dC, dV
    void compute_sensitivity(const std::vector<Eigen::VectorXd> &rho,
                             const std::vector<std::vector<Eigen::VectorXd>> &dH,
                             Eigen::VectorXd &dC, Eigen::VectorXd &dV);

    // simulation interface
    void simulation(const Eigen::MatrixXd &X, double &C, double &V, double &E, Eigen::VectorXd &dC,
                    Eigen::VectorXd &dV, Eigen::VectorXd &dE);

    // for point P with a global coordinate, get its responding element index
    // int getEleIDForPoint(const Eigen::RowVector3d &P) const;
    // update coordinates of points in triangle mesh (physical domain)
    void updateVInPD();
    void updateVInFineMesh();
    void output_fineCenter_stress(const std::vector<Eigen::VectorXd> &rho);
    void output_PD_stress();

    };

} // namespace SIM

#endif // MESH_HPP