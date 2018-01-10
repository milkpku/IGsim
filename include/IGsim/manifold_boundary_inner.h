/* This file is part of IGsim, a simple c++ simulation library.
 * 
 * Copyright (C) 2017 Like Ma <milkpku@gmail.com>
 * 
 * This Source Code Form is subject to the terms of the Mozilla Public License 
 * v. 2.0. If a copy of the MPL was not distributed with this file, You can 
 * obtain one at http://mozilla.org/MPL/2.0/.
 * This should *NOT* be contained in a IGSIM_*_H ifdef, since it may be defined
 * differently based on when it is included
 */

#ifndef IGSIM_MANIFOLD_BOUNDARY_INNER_H
#define IGSIM_MANIFOLD_BOUNDARY_INNER_H
#include "igsim_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sim
{
  /*  Construct boundary facet/edge and inner facet/edge for tetrahedron/triangular 
   *  meshes, and transfer attribute attached on tetrahedron/triangular to boundary
   *
   *  Inputs:
   *    T #T by {3|4} matrix of ordered tetrahedron/triangualr mesh
   *
   *  Outputs:
   *    F #F by {2|3} matrix of ordered face/edge, boundary of given T
   *    uF  #uF by {2|3} matrix of unordered face/edge, inner simplex of T
   */
  template<typename DerivedT, typename DerivedF, typename DeriveduF>
  IGSIM_INLINE void manifold_boundary_inner(
    const Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedF>& F, 
    Eigen::PlainObjectBase<DeriveduF>& uF);
  /*
   *  Inputs:
   *    TA #T by n matrix of information for tetrahedron/triangular mesh
   *
   *  Outputs:
   *    FA #F by n matrix of information for boundary face/line
   */
  template<
    typename DerivedT, typename DerivedTA, 
    typename DerivedF, typename DerivedFA>
  IGSIM_INLINE void manifold_boundary_inner(
    const Eigen::PlainObjectBase<DerivedT>& T,
    const Eigen::PlainObjectBase<DerivedTA>& TA,
    Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedFA>& FA);
}

#ifndef IGSTATIC_LIBRARY
#   include "manifold_boundary_inner.cpp"
#endif

#endif //!IGSIM_MANIFOLD_BOUNDARY_INNER_H
