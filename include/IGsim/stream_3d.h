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

#ifndef IGSIM_STREAM_3D_H
#define IGSIM_STREAM_3D_H
#include "igsim_inline.h"

#include <Eigen/Dense>


namespace sim
{
  /*  get streamlines of given vector field 
   *
   *  Inputs:
   *    V  #V by 3 matrix of coordinates of vertices
   *    T  #T by 4 matrix of indices of tetrahedron corners into vertices
   *    F  #T by 3 matrix of vector field in each tetrahedron
   *    dsep  double, separate 
   *
   *  Outputs:
   *    Ls #Ls by 3 matrix of coordinates of vertices on streamlines, is derived
   *       by concatenating all streamlines
   *    H  vector of accumulated length of streamlines in L, streamline i begins 
   *       at row H(i-1) and end before row H(i), 
   */
  template <
    typename DerivedV, typename DerivedT, typename DerivedF, 
    typename DerivedLs, typename DerivedH>
  IGSIM_INLINE void stream_3d(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    const Eigen::PlainObjectBase<DerivedF>& F,
    const double dsep,
    Eigen::PlainObjectBase<DerivedLs>& Ls,
    Eigen::PlainObjectBase<DerivedH>& H);

  /* get streamline of given vector field from a start point, given already
   * existing streamlines
   */
  template <
    typename DerivedV, typename DerivedT, typename DerivedF, 
    typename DerivedLs, typename DerivedL>
  IGSIM_INLINE void stream_3d(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    const Eigen::PlainObjectBase<DerivedF>& F,
    const Eigen::PlainObjectBase<DerivedP>& P,
    const Eigen::PlainObjectBase<DerivedLs>& Ls,
    const double dsep,
    Eigen::PlainObjectBase<DerivedL>& L);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "stream_3d.cpp"
#endif


#endif
