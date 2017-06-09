#include "python_shared.h"
#include <sstream>
#include <string>
#include <fstream>

extern void python_export_sim(py::module &);

PYBIND11_PLUGIN(pysim){
  py::module m("pysim", R"pysimdoc(
    Python wrapper for IGsim
    ------------------------
      average_onto_faces_mat
      boundary_facets
      manifold_harmonic_basis
      manifold_volume
      tet_topology
      writeVTK
  )pysimdoc");

  python_export_sim(m);

  return m.ptr();
}
