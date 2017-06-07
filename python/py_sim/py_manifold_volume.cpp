m.def("manifold_volume", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& T,
  Eigen::MatrixXd&  Vol
 )
{
  Eigen::VectorXd _Vol;
  sim::manifold_volume(V, T, _Vol);
  Vol = _Vol;
  return;
}, __doc_sim_manifold_volume,
py::arg("V"), py::arg("T"), py::arg("Vol")
);
