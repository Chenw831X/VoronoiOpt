# README

This project is the implementation of paper "Explicit Topology Optimization of Conforming Voronoi Foams".



The head files and library files are provided, along with an example of femur in paper. The library files were generated in **Linux_x86_64** platform.



### 1. Dependence

Our project is implemented in `C++` and built using `cmake`, having some dependences as follows:

* **TBB**
* **SuiteSparse (the CHOLMOD therein)**
* **Boost**
* **OpenMesh**
* **CGAL**
* **OpenMP**
* **geogram**
* **fmt**



### 2. Compile and Run

Run following commands in terminal

```shell
mkdir build && build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j16
./work
```



This will run the femur example set in `main.cpp`, and the following important settings in `main.cpp` are commented:

```
1. t, tmin, tmax: initial, minimum, maximum width of rods
2. scalarE: weight of shape energy
3. volfrac: volume fraction
```

The results are to be stored in the `output/femur` directory, including:

```
1. C.txt, E.txt, V.txt: compliance, shape energy and volume fraction of each iteration
2. model/model_i.txt: rod information of #iter-i, i.e. two endpoints and radius of rod
3. X/X_i.txt: varibles of #iter-i, i.e. seeds positions and radii
```



