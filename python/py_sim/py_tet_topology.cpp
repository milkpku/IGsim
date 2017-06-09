  
m.def("tet_topology", []
(
  const Eigen::MatrixXi& T,
  Eigen::MatrixXi& L
)
{
  return sim::tet_topology(T, L);
}, __doc_sim_tet_topology,
py::arg("T"), py::arg("L")
);

