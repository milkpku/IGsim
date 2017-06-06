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
