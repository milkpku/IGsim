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

#ifndef IGSIM_VOLUME_CLOSED_SURFACE_H
#define IGSIM_VOLUME_CLOSED_SURFACE_H
#include "igsim_inline.h"

#include <Eigen/Dense>

namespace sim
{
  /*  Calculate volume surrounded by closed surface manifold
   *  Inputs:
   *    V #V by 3 matrix of vertices coordinate
   *    F #F by 3 matrix of indices of vertices into V, representing manifold
   *      surface
   *
   *  Output:
   *    Vol scalar, volume surrounded by surface F, if F is not closed, the 
   *      the volume is sum of volume of tetrahedrons with geometry center as
   *      common vertices
   */
  template <typename DerivedV, typename DerivedF, typename ScalarVol>
  IGSIM_INLINE void volume_closed_surface(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    ScalarVol& Vol);
}

#ifndef IGSIM_STATIC_LIBRARY
  #include "volume_closed_surface.cpp"
#endif


#endif
