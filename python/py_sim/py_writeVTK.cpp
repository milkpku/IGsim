py::enum_<sim::VTK_TYPE>(m, "VTK_TYPE", py::arithmetic())
.value("VTK_SCALAR", sim::VTK_SCALAR)
.value("VTK_COLOR", sim::VTK_COLOR)
.value("VTK_VECTOR", sim::VTK_VECTOR)
.value("VTK_TEXTURE", sim::VTK_TEXTURE)
.export_values();

m.def("writeVTK", []
(
  const std::string& filename,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& T
)
{
  return sim::writeVTK(filename, V, T);
}, __doc_sim_writeVTK,
py::arg("filename"), py::arg("V"), py::arg("T")
);