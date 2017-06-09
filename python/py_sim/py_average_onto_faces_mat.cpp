  
m.def("average_onto_faces_mat", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::SparseMatrix<double>& Proj
)
{
  return sim::average_onto_faces_mat(V, F, Proj);
}, __doc_sim_average_onto_faces_mat,
py::arg("V"), py::arg("F"), py::arg("Proj")
);

