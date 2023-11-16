/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/4/20
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT__UTILS_H_
#define CVT__UTILS_H_

#include <igl/facet_components.h>
#include "common.h"
#include <igl/remove_unreferenced.h>
#include <fstream>
#include <iostream>
#include <deque>
#include <vector>
#include <string>


namespace cvt {
inline static void loadOpenMeshFromVF(TriMesh &mesh, const PntMat &V, const FaceMat &F){
    std::vector<OpenMesh::SmartVertexHandle> vhandles;
    for(int i = 0; i < V.rows(); i++){
        TriMesh::Point p(V(i, 0), V(i, 1), V(i, 2));
        vhandles.push_back(mesh.add_vertex(p));
    }
    for(int i = 0; i < F.rows(); i++){
        TriMesh::FaceHandle fh = mesh.add_face(vhandles[F(i, 0)], vhandles[F(i, 1)], vhandles[F(i, 2)]);
    }
}

inline static void dumpOpenMeshToVF(const TriMesh &mesh, PntMat &V, FaceMat &F){
    V.resize(mesh.n_vertices(), 3);
    F.resize(mesh.n_faces(), 3);
    for(auto vit : mesh.vertices()){
        V.row(vit.idx()) = Eigen::Vector3d(mesh.point(vit)[0], mesh.point(vit)[1], mesh.point(vit)[2]);
    }
    for(auto fit : mesh.faces()){
        auto vvec = fit.vertices().to_vector();
        F.row(fit.idx()) = Eigen::Vector3i(vvec[0].idx(), vvec[1].idx(), vvec[2].idx());
    }
}

inline static void separate(const PntMat &V, const FaceMat &F,
                            std::vector<PntMat> &VV, std::vector<FaceMat> &FF,
                            bool reindex=false){
    Eigen::MatrixXi C;
    igl::facet_components(F, C);
    std::vector<int> components_cnt(C.maxCoeff() + 1, 0);
    std::vector<int> components_row_pointers(C.maxCoeff() + 1, 0);
    VV.resize(C.maxCoeff() + 1);
    FF.resize(C.maxCoeff() + 1);
    for(int i = 0; i < C.rows(); i++){
        components_cnt[C(i)]++;
    }
    for (int i = 0; i < components_cnt.size(); ++i) {
        VV[i] = V;
        FF[i].resize(components_cnt[i], F.cols());
    }
    for (int i = 0; i < C.rows(); ++i) {
        int cid = C(i);
        int rid = components_row_pointers[cid];
        components_row_pointers[cid]++;
        FF[cid].row(rid) = F.row(i);
    }

    if(reindex){
        for(int i = 0; i < VV.size();++i){
            Eigen::MatrixXi J;
            PntMat VV_tmp = VV[i];
            FaceMat FF_tmp = FF[i];
            igl::remove_unreferenced(VV_tmp, FF_tmp, VV[i], FF[i], J);
        }
    }
}

void writeVTK(const std::string &path, const PntMat &V, const FaceMat &F, std::vector<double> facepatch);

void writeTetVTK(const std::string &path, const PntMat &V, const FaceMat &T, const std::vector<double> &cell_data = {}, const std::vector<double> &v_data = {});

VF readTetVTK(const std::string &path);

void writePntVTK(const std::string &path, const PntMat &V);

template<typename T>
void writeMatrix(const std::string &path, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M){
    std::ofstream ofs(path);
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            ofs << M(i, j) << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}

template<typename T>
void readMatrix(const std::string &path, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &M){
    std::vector<std::vector<T>> data;
    std::ifstream ifs(path);
    std::string line;
    int row = 0;
    int col = 0;
    while(std::getline(ifs, line)){
        std::stringstream ss(line);
        std::vector<T> row_data;
        col = 0;
        T v;
        while(ss >> v){
            col++;
            row_data.push_back(v);
        }
        data.push_back(row_data);
        row++;
    }

    M.resize(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            M(i, j) = data[i][j];
        }
    }
    ifs.close();
}
}

#endif //CVT__UTILS_H_
