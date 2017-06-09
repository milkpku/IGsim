  
m.def("manifold_harmonic_basis", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const int& k,
  Eigen::MatrixXcd& S,
  Eigen::MatrixXcd& U
)
{
  Eigen::VectorXcd _S;
  sim::manifold_harmonic_basis(V, F, k, _S, U);
  S = _S;
  return;
}, __doc_sim_manifold_harmonic_basis,
py::arg("V"), py::arg("F"), py::arg("k"), py::arg("S"), py::arg("U")
);

  
m.def("manifold_harmonic_basis", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const int& k,
  const double& sigma,
  Eigen::MatrixXcd& S,
  Eigen::MatrixXcd& U
)
{
  Eigen::VectorXcd _S;
  sim::manifold_harmonic_basis(V, F, k, sigma, _S, U);
  S = _S;
  return;
}, __doc_sim_manifold_harmonic_basis,
py::arg("V"), py::arg("F"), py::arg("k"), py::arg("sigma"), py::arg("S"), py::arg("U")
);

