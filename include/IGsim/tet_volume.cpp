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
#include "tet_volume.h"

#include <Eigen/Dense>


template <typename DerivedV, typename DerivedT, typename DerivedVol>
  IGSIM_INLINE void sim::tet_volume(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    Eigen::PlainObjectBase<DerivedVol>& Vol)
{
  assert(V.cols() == 3 && "V.cols() should be 3");
  assert(T.cols() == 4 && "T should be tetrahedron mesh");

  Vol.resize(T.rows());

  for (int i = 0; i < T.rows(); i++)
  {
    auto t = T.row(i);
    
    Eigen::MatrixXd Dm(3, 3);
    Dm.row(0) << V.row(t(0)) - V.row(t(3));
    Dm.row(1) << V.row(t(1)) - V.row(t(3));
    Dm.row(2) << V.row(t(2)) - V.row(t(3));

    Vol(i) = Dm.determinant();
  }
  Vol /= 6.0;
}

