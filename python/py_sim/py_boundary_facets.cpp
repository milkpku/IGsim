m.def("boundary_facets", []
(
  const Eigen::MatrixXi& T,
  Eigen::MatrixXi&       F
)
{
  return sim::boundary_facets(T, F);
}, __doc_sim_boundary_facets,
py::arg("T"), py::arg("F")
);

m.def("boundary_facets", []
(
  const Eigen::MatrixXi& T,
  const Eigen::MatrixXd& TA,
  Eigen::MatrixXi&       F,
  Eigen::MatrixXd&       FA
)
{
  return sim::boundary_facets(T, TA, F, FA);
}, __doc_sim_boundary_facets,
py::arg("T"), py::arg("TA"), py::arg("F"), py::arg("FA")
);
