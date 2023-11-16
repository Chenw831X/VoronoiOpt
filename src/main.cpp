#include <iostream>
#include <chrono>
#include <Eigen/Eigen>
#include <spdlog/spdlog.h>
#include <igl/combine.h>

#include "CVT/CVT.h"
#include "CVT/utils.h"

#include "Type.hpp"
#include "Utils.hpp"
#include "Mesh.hpp"

#include <mma/MMASolver.h>

std::shared_ptr<SIM::PhysicalDomain> physicalDomain;
std::shared_ptr<cvt::CVTAlgAdvanced> cvt_alg;
std::shared_ptr<SIM::Mesh> mesh;
std::shared_ptr<MMASolver> mma;
int nCell, nVar;

std::map<std::string, double> timer;

/* config */
int example;

// triangle mesh file path of physical domain
std::string workspace;
std::string PD_FilePath;
std::string cvt_mesh;
std::string seeds_path;
std::string macro_prefix;
int macro_size;
std::string tet_prefix;

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
/* config */

void init_config() {
    example = 2; // femur

    // Material Setup
    YM1 = 1.0e5;
    YM0 = 1.0e-5;
    PR = 0.3;
    penaltyYM = 1;
    penaltyDBC = 1e10;

    workspace = "../output/femur/";
    std::string model = "femur";

    cvt_mesh = workspace + "input/femur.obj";
    PD_FilePath = workspace + "input/femur.obj";
    macro_prefix = workspace + "macro/cell";
    macro_size = 640;

    // Boundary Condition Setup
    // NBC (downward)
    double NBC_val = 10.0;
    Eigen::Matrix<double, 2, 3> tempNBCRelBBox;
    Eigen::Vector3d tempNBCVal;
    // NBC-1
    tempNBCRelBBox.row(0) = Eigen::RowVector3d(0.1550, 0.9807, 0.5673);
    tempNBCRelBBox.row(1) = Eigen::RowVector3d(0.3823, 1.01, 0.8343);
    tempNBCVal << 0.0, -1.0, 0.0;
    tempNBCVal *= NBC_val;
    NBCRelBBox.emplace_back(tempNBCRelBBox);
    NBCVal.emplace_back(tempNBCVal);

    // NBC-2
    tempNBCRelBBox.row(0) = Eigen::RowVector3d(0.7767, 0.7206, 0.1525);
    tempNBCRelBBox.row(1) = Eigen::RowVector3d(0.9175, 0.7682, 0.3490);
    tempNBCVal << 0.0, 1.0, 0.0;
    tempNBCVal *= NBC_val;
    NBCRelBBox.emplace_back(tempNBCRelBBox);
    NBCVal.emplace_back(tempNBCVal);

    Eigen::Matrix<double, 2, 3> tempDBCRelBBox;
    tempDBCRelBBox.row(0) = Eigen::RowVector3d(-0.1, -0.1, -0.1);
    tempDBCRelBBox.row(1) = Eigen::RowVector3d(1.1, 0.001, 1.1);
    DBCRelBBox.emplace_back(tempDBCRelBBox);
    DBCVal.emplace_back(Eigen::Vector3d(0.0, 0.0, 0.0));

    // t = 0.06; tmin = 0.042; tmax = 0.078; // 100 seeds
    t = 0.026; tmin = 0.013; tmax = 0.039; // 500 seeds
    // t = 0.018; tmin = 0.009; tmax = 0.027; // 5000 seeds
    // scalarE = 0.01;
    scalarE = 0.5;
    // scalarE = 0.1;
    volfrac = 0.3;

    seeds_path = workspace + "seeds.txt";
    tet_prefix = workspace + "tets/tet";

}

void display_model(const std::string &path){
    auto edge_vec = cvt_alg->voronoi_mesh.voro_edge_vectors();
    // auto edge_vec = cvt_alg->voronoi_mesh.voro_detail_edge_vectors();
    int sz = (int)edge_vec.size();
    Eigen::MatrixXd model(sz, 7);
    int cnt = 0;

    for(const auto &edge : edge_vec){
        model(cnt, 0) = std::get<3>(edge);
        model(cnt, 1) = std::get<2>(edge).first(0);
        model(cnt, 2) = std::get<2>(edge).first(1);
        model(cnt, 3) = std::get<2>(edge).first(2);
        model(cnt, 4) = std::get<2>(edge).second(0);
        model(cnt, 5) = std::get<2>(edge).second(1);
        model(cnt, 6) = std::get<2>(edge).second(2);
        ++cnt;
    }
    assert(cnt == sz);
    SIM::Utils::writeMatrixXd(path, model);
}

void initialize() {
    auto preprocess_start = std::chrono::steady_clock::now();

    /* initialize triangle mesh */
    physicalDomain = std::make_shared<SIM::PhysicalDomain>(PD_FilePath,
                                                           NBCRelBBox, NBCVal,
                                                           DBCRelBBox, DBCVal);
#ifdef FCM_NBC
    physicalDomain->writeNBC();
#endif

#ifdef FCM_DBC
    physicalDomain->writeDBC_obj();
#endif

    /* initialize cvt */
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(cvt_mesh, V, F);
    cvt_alg = std::make_shared<cvt::CVTAlgAdvanced>(V, F);
    // cvt_alg = std::make_shared<cvt::CVTAlg>(V, F);

    if (example == 1){
        cvt_alg->background_mesh(
                cvt::seed::FileSeedLoader(seeds_path),
                cvt::tet::TetMeshLoader(tet_prefix),
                t          // radius
        );
    }
    else{
        cvt_alg->background_mesh(macro_prefix,
                                 macro_size,
                                 cvt::tet::TetMeshLoader(tet_prefix),
                                 t);
    }
    Eigen::Vector3d len = cvt_alg->box.sizes();
    cvt_alg->box.min() -= len * 0.2;
    cvt_alg->box.max() += len * 0.2;
    std::cout << "------ box ------" << std::endl;
    std::cout << cvt_alg->box.min().transpose() << std::endl;
    std::cout << cvt_alg->box.max().transpose() << std::endl;

    // compute polygon, micV, micT
    int nEle = cvt_alg->tet_mesh.size();
    assert(cvt_alg->bg_voronoi_mesh.cells_connectivities.size() == nEle);
    std::vector<std::vector<Eigen::MatrixXd>> polygon(nEle);
    std::vector<Eigen::MatrixXd> micV(nEle);
    std::vector<Eigen::MatrixXi> micT(nEle);
    for(int macI=0; macI<nEle; ++macI){
        micV[macI] = cvt_alg->tet_mesh[macI].TV;
        micT[macI] = cvt_alg->tet_mesh[macI].TT;

        const auto &cell = cvt_alg->bg_voronoi_mesh.cells_connectivities[macI].cell;
        int nFace = cell.patches.size();
        // std::cout << "nFace: " << nFace << std::endl;
        polygon[macI].resize(nFace);
        for(int fI=0; fI<nFace; ++fI){
            int nSeg = cell.patches[fI].edges_in_order.size();
            polygon[macI][fI].resize(nSeg, 6);
            for(int sI=0; sI<nSeg; ++sI){
                const auto &edge = cell.patches[fI].edges_in_order[sI];
                const int &front_ = edge.first;
                const int &end_ = edge.second;
                auto pf = cell.mesh.point(cell.mesh.vertex_handle(front_));
                polygon[macI][fI](sI, 0) = pf[0];
                polygon[macI][fI](sI, 1) = pf[1];
                polygon[macI][fI](sI, 2) = pf[2];
                auto pb = cell.mesh.point(cell.mesh.vertex_handle(end_));
                polygon[macI][fI](sI, 3) = pb[0];
                polygon[macI][fI](sI, 4) = pb[1];
                polygon[macI][fI](sI, 5) = pb[2];
            }
        }
    }

    // use shell
    if (example == 2) {
        cvt_alg->set_shell(true, 0.03);
    }
    else if (example == 4){
        cvt_alg->set_shell(true, 0.0015);
    }

    if(example == 1){
        Eigen::MatrixXd X = cvt_alg->get_vars();
        nCell = X.rows();
        cvt_alg->step(X); // need to step in the first
        nVar = nCell * 4;
    }
    else{
        Eigen::MatrixXd X0;
        cvt::readMatrix(seeds_path, X0);
        nCell = X0.rows();
        assert(X0.cols() == 3);
        X0.conservativeResize(nCell, 4);
        X0.col(3).setConstant(t);
        cvt_alg->step(X0);
        nVar = nCell * 4;
    }
    spdlog::info("number of voronoi cell: {}", nCell);


    /* initialize CBN mesh */
    mesh = std::make_shared<SIM::Mesh>(YM1, YM0, PR, penaltyYM,
                                       penaltyDBC, nCell, nVar, physicalDomain, cvt_alg,
                                       polygon, micV, micT);

    std::cout << "step: " << mesh->step << std::endl;

    mesh->computePhi();
    mesh->preprocess();
    mesh->prepareMicroIntegration();

    /* initialize MMA */
    int nCon = 1;
    mma = std::make_shared<MMASolver>(nVar, nCon);

    std::chrono::duration<double> preprocess_(std::chrono::steady_clock::now() - preprocess_start);
    timer["preprocess"] += preprocess_.count();
    spdlog::info("preprocess time: {:.6f}", timer["preprocess"]);

}


void optimization() {
    int dgt0 = 5;
    int loop = 0, maxloop = 100;
    std::vector<double> OBJ;
    double objVr5 = 1.0;

    Eigen::MatrixXd X = cvt_alg->get_vars();
    // set min, max of variable
    Eigen::VectorXd Xmin(nVar);
    Eigen::VectorXd Xmax(nVar);
    Eigen::Matrix<double, 2, 3> meshBBox;
    meshBBox.row(0) = mesh->node.colwise().minCoeff();
    meshBBox.row(1) = mesh->node.colwise().maxCoeff();
    std::cout << "meshBBox: " << std::endl << meshBBox << std::endl;
    Xmin(Eigen::seq(0, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, meshBBox(0, 0));
    Xmin(Eigen::seq(1, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, meshBBox(0, 1));
    Xmin(Eigen::seq(2, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, meshBBox(0, 2));
    Xmin(Eigen::seq(3, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, tmin);
    Xmax(Eigen::seq(0, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, meshBBox(1, 0));
    Xmax(Eigen::seq(1, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, meshBBox(1, 1));
    Xmax(Eigen::seq(2, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, meshBBox(1, 2));
    Xmax(Eigen::seq(3, Eigen::last, 4)) = Eigen::VectorXd::Constant(nCell, tmax);

    // for output
    display_model(workspace + "model/model_0.txt");
    Eigen::VectorXd C_out, V_out, E_out;

    // while (loop < maxloop && objVr5 > 1.0e-3) {
    while (loop < maxloop) {
        auto iter_start = std::chrono::steady_clock::now();

        ++loop;
        double C, V, E;
        Eigen::VectorXd dC;
        Eigen::VectorXd dV;
        Eigen::VectorXd dE;
        mesh->simulation(X, C, V, E, dC, dV, dE);
        spdlog::debug("dC.mean = {}", dC.mean());
        spdlog::debug("dV.mean = {}", dV.mean());
        spdlog::debug("dE.mean = {}", dE.mean());

        int dgt_c = dgt0 - int(log10(dC.cwiseAbs().maxCoeff()));
        int dgt_v = dgt0 - int(log10(dV.cwiseAbs().maxCoeff()));
        int dgt_e = dgt0 - int(log10(dE.cwiseAbs().maxCoeff()));
        double tmp_c = pow(10.0, (double)dgt_c);
        double tmp_v = pow(10.0, (double)dgt_v);
        double tmp_e = pow(10.0, (double)dgt_e);
        dC = (dC * tmp_c).array().round() / tmp_c;
        dV = (dV * tmp_v).array().round() / tmp_v;
        dE = (dE * tmp_e).array().round() / tmp_e;

        // normalization
        dC /= dC.maxCoeff();
        dV /= dV.maxCoeff();
        dE /= dE.maxCoeff();

        Eigen::VectorXd X_tmp = X.transpose().reshaped(nVar, 1);
        double f0val = C * (1.0 - scalarE) + E * scalarE;
        Eigen::VectorXd df0dx = -dC * (1.0 - scalarE) + dE * scalarE;
        double fval = V / volfrac - 1;
        Eigen::VectorXd dfdx = dV / volfrac;

        spdlog::info("mma update");
        mma->Update(X_tmp.data(), df0dx.data(), &fval, dfdx.data(), Xmin.data(), Xmax.data());
        X = X_tmp.reshaped(4, nCell).transpose();
        spdlog::debug("(updated) X.mean = {}", X.mean());
        spdlog::info("step: update design variable");

        auto bool_start = std::chrono::steady_clock::now();
        cvt_alg->step(X); // update design variable
        std::chrono::duration<double> bool_(std::chrono::steady_clock::now() - bool_start);
        timer["bool"] += bool_.count();

        // compute f0val
        OBJ.emplace_back(f0val);
        if(loop >= 5 && (V - volfrac) / volfrac < 1.0e-3){
            double mean_ = 0.0;
            for(int i_=loop-5; i_<loop; ++i_){
                mean_ += OBJ[i_];
            }
            mean_ /= 5.0;
            double max_ = 0.0;
            for(int i_=loop-5; i_<loop; ++i_){
                max_ = std::max(max_, abs(OBJ[i_] - mean_));
            }
            objVr5 = abs(max_ / mean_);
        }

        // output
        C_out.conservativeResize(loop);
        C_out(loop - 1) = C;
        V_out.conservativeResize(loop);
        V_out(loop - 1) = V;
        E_out.conservativeResize(loop);
        E_out(loop - 1) = E;
        physicalDomain->writeOBJ(loop);
        SIM::Utils::writeMatrixXd(workspace + "fine_V/fine_V_" + std::to_string(loop) + ".txt", mesh->fine_V_deformed);
        SIM::Utils::writeMatrixXd(workspace + "X/X_" + std::to_string(loop) + ".txt", X);
        display_model(workspace + "model/model_" + std::to_string(loop) + ".txt");
        spdlog::critical("Optimization iter# {}, C = {:.6}, E = {:.6f}, V = {:.3f}, ch = {:.4f}", loop, C, E, V, objVr5);

        std::chrono::duration<double> iter_(std::chrono::steady_clock::now() - iter_start);
        timer["iter"] += iter_.count();

        spdlog::info("until now: total simulation time: {:.6f}", timer["sim"]);
        spdlog::info("until now: total dH time: {:.6f}", timer["dH"]);
        spdlog::info("until now: total bool time: {:.6f}", timer["bool"]);
        spdlog::info("until now: total rho time: {:.6f}", timer["rho"]);
        spdlog::info("until now: total iter time: {:.6f}", timer["iter"]);
    }

    SIM::Utils::writeMatrixXd(workspace + "C.txt", C_out);
    SIM::Utils::writeMatrixXd(workspace + "V.txt", V_out);
    SIM::Utils::writeMatrixXd(workspace + "E.txt", E_out);

}

int main(int argc, char **argv) {
    auto totalStart = std::chrono::steady_clock::now();

    // initialize
    // 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, 6=off
    spdlog::set_level(spdlog::level::trace); // print all logs
#ifdef USE_TBB
    int num_threads = oneapi::tbb::info::default_concurrency();
    spdlog::info("TBB: number of threads: {}", num_threads);
#endif
    Eigen::setNbThreads(8);
    spdlog::info("OMP/Eigen: number of threads: {}", Eigen::nbThreads());

    // timer
    timer["preprocess"] = 0.0;
    timer["iter"] = 0.0;
    timer["sim"] = 0.0;
    timer["dH"] = 0.0;
    timer["bool"] = 0.0;
    timer["rho"] = 0.0;

    init_config();
    initialize();
    optimization();

    std::chrono::duration<double> total_(std::chrono::steady_clock::now() - totalStart);
    spdlog::info("program finish, total time cost: {:.6f} s", total_.count());
    return 0;
}