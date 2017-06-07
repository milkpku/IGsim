#include <Eigen/Dense>

#include "python_shared.h"

#include <IGsim/boundary_facets.h>
#include <IGsim/manifold_volume.h>
#include <IGsim/writeVTK.h>

void python_export_sim(py::module &m)
{
  #include "py_sim/py_boundary_facets.cpp"
  #include "py_sim/py_manifold_volume.cpp"
  #include "py_sim/py_writeVTK.cpp"
}
