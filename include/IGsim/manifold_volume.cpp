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

#include "manifold_volume.h"

#include <igl/doublearea.h>

template <typename DerivedV, typename DerivedF, typename DerivedVol>
IGSIM_INLINE void sim::manifold_volume(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedF>& F,
  Eigen::PlainObjectBase<DerivedVol>& Vol)
{
  if (F.cols() == 3)
  {
    igl::doublearea(V, F, Vol);
    Vol /= 2;
  }
  else if (F.cols() == 4)
  {
    assert(V.cols() == 3 && "only support 3D vertices");
    Vol = DerivedVol::Zero(F.rows());
	  for (int i = 0; i < F.rows(); i++)
	  {
		  Eigen::Matrix3d tet;
		  tet << V.row(F(i, 0)), V.row(F(i, 1)), V.row(F(i, 2));
		  tet.rowwise() -= V.row(F(i, 3));
		  Vol(i) = tet.determinant() / 6;
	  }
  }
  else
  {
    assert(false && "only 2D and 3D manifold supported");
  }
}
