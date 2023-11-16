/**
 * ------------------------------------
 * @author: Weipeng Kong
 * @date: 2022/5/6
 * @email: yjxkwp\@foxmail.com
 * @description: 
 * ------------------------------------
**/

#ifndef CVT_INCLUDE_CVT_CVTSEEDS_H_
#define CVT_INCLUDE_CVT_CVTSEEDS_H_

#include <Eigen/Eigen>
#include <vector>
#include <boost/filesystem.hpp>
#include <memory>

namespace cvt {
namespace seed {

class SeedGenerator {
public:
    virtual std::vector<Eigen::Vector3d> operator()(Eigen::AlignedBox<double, 3> &bbox,
                                                    const Eigen::MatrixXd &V,
                                                    const Eigen::MatrixXi &F) const = 0;
};

class RandomSeedGenerator : public SeedGenerator {
    int num_seeds;
public:
    explicit RandomSeedGenerator(int num_seeds);

    std::vector<Eigen::Vector3d> operator()(Eigen::AlignedBox<double, 3> &bbox,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F) const;
};

class ImportanceSamplingSeedGenerator : public SeedGenerator {
    int uniform_sample_n;
    int corner_sample_n;
    int inmesh_sample_n;
    double inmesh_sample_iso;
    double outer_scale_radio;
    double inner_scale_radio;

public:
    ImportanceSamplingSeedGenerator(int uniform_sample_n,
                                    int corner_sample_n,
                                    int inmesh_sample_n,
                                    double inmesh_sample_iso,
                                    double outer_scale_radio,
                                    double inner_scale_radio);

    std::vector<Eigen::Vector3d> operator()(Eigen::AlignedBox<double, 3> &bbox,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F) const;
};

class FileSeedLoader : public SeedGenerator {
    boost::filesystem::path file_path;
public:
    FileSeedLoader(const boost::filesystem::path &file_path);
    std::vector<Eigen::Vector3d> operator()(Eigen::AlignedBox<double, 3> &bbox,
                                            const Eigen::MatrixXd &V,
                                            const Eigen::MatrixXi &F) const;
};

void SaveSeeds(const boost::filesystem::path &file_name, const Eigen::AlignedBox<double, 3> &bbox, const std::vector<Eigen::Vector3d> &seeds);

};
}

#endif //CVT_INCLUDE_CVT_CVTSEEDS_H_
