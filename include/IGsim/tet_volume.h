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

#ifndef IGSIM_TET_VOLUME_H
#define IGSIM_TET_VOLUME_H
#include "igsim_inline.h"

#include <Eigen/Dense>

namespace sim
{
  /*  Compute tetrahedron volume for given tetmesh, each tetrahedron is
   *  [t0, t1, t2, t3] and its boundary towards outside is 
   *  [t0, t1, t2] - [t1, t2, t3] + [t2, t3, t0] - [t3, t0, t1]
   *
   *  Inputs:
   *    V #V by 3 matrix of coordinate of vertices
   *    T #T by 4 matrix of tetrahedron mesh
   * 
   *  Output:
   *    Vol #T vector of tetrahedron volume, with direction
   */
  template <typename DerivedV, typename DerivedT, typename DerivedVol>
    IGSIM_INLINE void tet_volume(
      const Eigen::PlainObjectBase<DerivedV>& V,
      const Eigen::PlainObjectBase<DerivedT>& T,
      Eigen::PlainObjectBase<DerivedVol>& Vol);
}


#ifndef IGSIM_STATIC_LIBRARY
#   include "tet_volume.cpp"
#endif

#endif
