#include <Eigen/Dense>

#include "python_shared.h"

#include <IGsim/boundary_facets.h>

void python_export_sim(py::module &m)
{
  #include "py_sim/py_boundary_facets.cpp"
}
