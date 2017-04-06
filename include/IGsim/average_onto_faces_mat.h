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

#ifndef   IGSIM_AVERAGE_ONTO_FACES_MAT_H
#define   IGSIM_AVERAGE_ONTO_FACES_MAT_H
#include "igsim_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sim
{
  /*  Construct the average matrix from vertices to faces, which turn nodal
   *  scalar to facial scalar
   *
   *  Inputs:
   *    V #V by x matrix of vertex coordinates
   *    F #F by x matrix, indices of simplex vertex into V
   *
   *  Outputs:
   *    Proj  #F by #V sparse matrix, SF = Proj * SV, where SF and SV are scalar
   *      represented on faces and vertices respectively
   */
  template <typename DerivedV, typename DerivedF, typename Scalar>
  IGSIM_INLINE void average_onto_faces_mat(
      const Eigen::PlainObjectBase<DerivedV>& V,
      const Eigen::PlainObjectBase<DerivedF>& F, 
      Eigen::SparseMatrix<Scalar>& Proj);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "average_onto_faces_mat.cpp"
#endif

#endif //!IGSIM_AVERAGE_ONTO_FACES_MAT_H

