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
#include "tet_topology.h"

#include <Eigen/Dense>

#include <igl/sortrows.h>

template <typename DerivedT, typename DerivedL>
IGSIM_INLINE void sim::tet_topology(
  const Eigen::PlainObjectBase<DerivedT>& T,
  Eigen::PlainObjectBase<DerivedL>& L)
{
  DerivedT TF( 4 * T.rows(), 5);

  /* construct undirected surface and its associate tet info */
  for(int i = 0; i < T.rows(); i++)
  {
    Eigen::VectorXi  id;
  DerivedT unsort_t, t;
  unsort_t = T.row(i).transpose();
    igl::sortrows(unsort_t, true, t, id);

    TF.row(4 * i    ) << t(1), t(2), t(3), i, id(0);
    TF.row(4 * i + 1) << t(0), t(2), t(3), i, id(1);
    TF.row(4 * i + 2) << t(0), t(1), t(3), i, id(2);
    TF.row(4 * i + 3) << t(0), t(1), t(2), i, id(3);
  }

  /* sort tet surface to get adjacent tet together */
  DerivedT sort_TF(TF.rows(), TF.cols());
  Eigen::VectorXi id_TF;
  igl::sortrows(TF, true, sort_TF, id_TF);

  /* construct topology matrix L */
  L = - DerivedL::Ones(T.rows(), T.cols());
  DerivedT surf = sort_TF.leftCols(3);
  DerivedT info = sort_TF.rightCols(2);
  int count = 1;
  while (count < sort_TF.rows())
  {
    if ((surf.row(count).array() == surf.row(count-1).array()).all())
    {
      L(info(count, 0), info(count, 1)) = info(count - 1, 0);
      L(info(count - 1, 0), info(count - 1, 1)) = info(count, 0);
    count += 2;
    }
  else {
    count++;
  }
  }
}
