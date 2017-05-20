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

#ifndef IGSIM_TET_TOPOLOGY_H
#define IGSIM_TET_TOPOLOGY_H
#include "igsim_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sim
{
  /*  Construct the topology of tetrahedrons
   *
   *  Inputs:
   *    T   #T by 4 matrix of indices of vertices of tetrahedron corners, 
   *        ordered by [v1, v2, v3, v4], and the boudanry aim towards outside.
   *
   *  Outputs:
   *    L   #T by 4 matrix of indices of neighbor tetrahedrons correspond to 
   *        vertices [v1, v2, v3, v4], or -1 when no neighbor tetrahedron   
   */
  template <typename DerivedT, typename DerivedL>
  IGSIM_INLINE void tet_topology(
    const Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedL>& L);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "tet_topology.cpp"
#endif


#endif
