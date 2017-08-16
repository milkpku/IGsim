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

#include "volume_closed_surface.h"

#include <Eigen/Dense>

template <typename DerivedV, typename DerivedF, typename ScalarVol>
IGSIM_INLINE void sim::volume_closed_surface(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  ScalarVol& Vol)
{
  assert(V.cols() == 3 && "should be 3D vertices");
  assert(F.cols() == 3 && "should be triangular mesh");

  /* gravity/geometry center */
  const Eigen::RowVector3d o = V.colwise().mean();
  
  Vol = 0;
  for (int i = 0; i < F.rows(); ++i)
  {
    Eigen::RowVector3d a = V.row(F(i, 0));
    Eigen::RowVector3d b = V.row(F(i, 1));
    Eigen::RowVector3d c = V.row(F(i, 2));

    Vol += (a - o).dot((b - o).cross(c - o));
  }

  Vol /= 6.;
}
