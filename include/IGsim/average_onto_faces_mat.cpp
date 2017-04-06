/* This file is part of IGsim, a simple c++ simulation library.
 * 
 * Copyright (C) 2016 Like Ma <milkpku@gmail.com>
 * 
 * This Source Code Form is subject to the terms of the Mozilla Public License 
 * v. 2.0. If a copy of the MPL was not distributed with this file, You can 
 * obtain one at http://mozilla.org/MPL/2.0/.
 * This should *NOT* be contained in a IGSIM_*_H ifdef, since it may be defined
 * differently based on when it is included
 */

#include "average_onto_faces_mat.h"


template <typename DerivedV, typename DerivedF, typename Scalar>
IGSIM_INLINE void sim::average_onto_faces_mat(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F, 
  Eigen::SparseMatrix<Scalar>& Proj)
{
  Proj.setZeros();
  Proj.resize(F.rows(), V.rows());

  typedef Eigen::Triplet<Scalar> T;
  std::vector<T> proj_coeff;
  proj_coeff.clear();
  proj_coeff.reserve(F.size());

  const Scalar avg = 1.0/F.cols();

  for(int j = 0; j < F.cols(); j++)
    for(int i = 0; i < F.rows(); i++)
      proj_coeff.push_back(T(i, F(i, j), avg));

  Proj.setfromTriplets(proj_coeff.begin(), proj_coeff.end());
}
