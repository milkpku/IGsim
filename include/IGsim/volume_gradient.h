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

#ifndef IGSIM_VOLUME_GRADIENT_H
#define IGSIM_VOLUME_GRADIENT_H
#include "igsim_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sim
{
  /*  Calculate volume gradient of volume surrounded by closed surface manifold, and its volume 
   *  gradient, and its Hessian
   *
   *  Inputs:
   *    V #V by 3 matrix of vertices coordinate
   *    F #F by 3 matrix of indices of vertices into V, representing manifold
   *      surface
   *
   *  Outputs:
   *    G #V by 3, volume gradient on each vertex
   */
  template <typename DerivedV, typename DerivedF, typename DerivedG>
  IGSIM_INLINE void volume_gradient(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedG>& G);
  /*  Outputs:
   *    H (3 * #V) by (3 * #V) hessian of volume
   */
  template <typename DerivedV, typename DerivedF, typename ScalarH>
  IGSIM_INLINE void volume_gradient(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::SparseMatrix<ScalarH>& H);
}

#ifndef IGSIM_STATIC_LIBRARY
  #include "volume_gradient.cpp"
#endif

#endif
