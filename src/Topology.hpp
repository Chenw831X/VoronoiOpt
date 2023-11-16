//
// Created by Wei Chen on 3/4/22
//

#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include <Eigen/Eigen>

namespace SIM{

// base class of topology components
class Topology{
public:
    Eigen::MatrixXd Coords;
    int nNodeOneDim;

public:
    int nNode;
    Eigen::VectorXi vIndex; // index of nodes
    Eigen::VectorXi dofIndex; // index of dofs
    Eigen::MatrixXd V0; // coordinates of nodes

public:
    Topology(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim);

public:
    virtual int assignVIndex(int idx); // assign dof
};

// class Topology Node and Its Geometry Information
class Node : public Topology{
    typedef Topology Base;

public:
    Node(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim);
};

// class Topology Edge and Its Geometry Information
// A Edge has a Vertex distribution as follows:
//  0        1
//  o--------o
// 01 is x-direction, y-direction, or z-direction
class Edge : public Topology{
    typedef Topology Base;

public:
    Edge(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim);
};

// class Topology Face and Its Geometry Information
// A face has a Vertex distribution as follows:
//  3        2
//  o--------o
//  |        |
//  |        |
//  |        |
//  o--------o
//  0        1
// 01 is the direction 1, 03 is the direction 2
// only 3 cases:
//  xy face: 01 is x-direction, 03 is y-direction
//  yz face: 01 is y-direction, 03 is z-direction
//  zx face: 01 is z-direction, 03 is x-direction
class Face : public Topology{
    typedef Topology Base;

public:
    Face(const Eigen::MatrixXd &p_Coords, int p_nNodeOneDim);
};

} // namespace SIM

#endif // TOPOLOGY_HPP