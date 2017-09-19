# IGsim - A tiny C++ FEM simulation library

IGsim is a tiny C++ FEM simulation library. It has implemented functions which calculate deformation energy, force and sitffness matrix of neo-hookean material and some other tools which are handy for 3D FEM based simulation.

It is a **header-only library**, mainly mimic the design philosophy of [libigl](https://github.com/libigl/libigl/). You do not need to compile anything to use, just include IGsim headers (e.g. `#include <IGsim/elastic_neohookean.h>`) and run. Most are tailored to operate on a generic triangle mesh stored in an n-by-3 matrix of vertex positions `V`, an m-by-3 matrix of triangle indices `F` or an m-by-4 matrix of tetrahedron indices `T`.

IGsim uses the [Eigen](http://eigen.tuxfamily.org) library heavily, mainly as input/output data container and matrix algebra interface. The representation of meshes is coordinated with libigl, and you can refer to libigl's [tutorial](http://libigl.github.io/libigl/tutorial/tutorial.html) for more information.
