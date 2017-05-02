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

#ifndef IGSIM_MANIFOLD_VOLUME_H
#define IGSIM_MANIFOLD_VOLUME_H
#include "igsim_inline.h"

#include <Eigen/Dense>

namespace sim
{
  /*  Calculate the manifold volume for each simplex (facet area for each 
   *  triangle in 2D manifold and volume for each tetrahedron in 3D manifold)
   *
   *  Inputs:
   *    V #V by dim matrix of vertex coordinates
   *    F #F by {3|4} matrix of indices of simplex vertex into V
   *
   *  Output:
   *    Vol #F by 1 matrix of manifold volume
   */
  template <typename DerivedV, typename DerivedF, typename DerivedVol>
  IGSIM_INLINE void manifold_volume(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    Eigen::PlainObjectBase<DerivedVol>& Vol);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "manifold_volume.cpp"
#endif

#endif //!IGSIM_MANIFOLD_VOLUME_H

