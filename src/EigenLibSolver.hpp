//
//  Created by Wei Chen on 3/31/22
//

#ifndef EigenLibSolver_hpp
#define EigenLibSolver_hpp

#include "LinSysSolver.hpp"

#include <Eigen/Eigen>

#include <vector>
#include <set>

namespace SIM {

template <typename vectorTypeI, typename vectorTypeS>
class EigenLibSolver : public LinSysSolver<vectorTypeI, vectorTypeS> {
    typedef LinSysSolver<vectorTypeI, vectorTypeS> Base;

protected:
    Eigen::SparseMatrix<double> coefMtr;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> simplicialLDLT;

public:
    void set_pattern(const std::vector<std::set<int>>& vNeighbor);
    void set_pattern(const Eigen::SparseMatrix<double>& mtr); //NOTE: mtr must be SPD

    void analyze_pattern(void);

    bool factorize(void);

    void solve(Eigen::VectorXd& rhs,
        Eigen::VectorXd& result);

    double coeffMtr(int rowI, int colI) const;

    void setZero(void);

    virtual void setCoeff(int rowI, int colI, double val);

    virtual void addCoeff(int rowI, int colI, double val);

    virtual void setUnit_row(int rowI);

    virtual void setUnit_col(int colI, const std::set<int>& rowVIs);
};

} // namespace SIM

#endif /* EigenLibSolver_hpp */
