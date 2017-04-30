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

#include "igsim_inline.h"
#include "stream_3d.h"

#include <vector>
#include <queue>

template <
  typename DerivedV, typename DerivedT, typename DerivedF, 
  typename DerivedLs, typename DerivedH>
IGSIM_INLINE void sim::stream_3d(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedT>& T,
  const Eigen::PlainObjectBase<DerivedF>& F,
  const double dsep,
  Eigen::PlainObjectBase<DerivedLs>& Ls,
  Eigen::PlainObjectBase<DerivedH>& H)
{
  /*TODO stream algorithm*/
  
}

template <
  typename DerivedV, typename DerivedT, typename DerivedF, 
  typename DerivedLs, typename DerivedL>
IGSIM_INLINE void sim::stream_3d(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedT>& T,
  const Eigen::PlainObjectBase<DerivedF>& F,
  const Eigen::PlainObjectBase<DerivedP>& P,
  const Eigen::PlainObjectBase<DerivedLs>& Ls,
  const double dsep,
  Eigen::PlainObjectBase<DerivedL>& L)
{
  /*TODO stream algorithm*/
}
