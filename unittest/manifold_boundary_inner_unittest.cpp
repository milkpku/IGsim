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

#include "../include/IGsim/manifold_boundary_inner.h"
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <igl/readOBJ.h>

#include "unittest_defines.h"

typedef Eigen::VectorXi iVec;
typedef Eigen::VectorXd dVec;
typedef Eigen::MatrixXd dMat;
typedef Eigen::MatrixXi iMat;

TEST(manifold_boundary_inner, facet_num)
{
  /* tet cube 
   *   4 - - - 7
   *   |\      |\
   *   | 5 - - - 6
   *   0 | - - 3 |
   *    \|      \|
   *     1 - - - 2
   * */
  iMat T(5, 4);
  T << 0, 1, 4, 3,
       2, 3, 6, 1,
       4, 5, 6, 1,
       4, 7, 6, 3,
       1, 3, 6, 4;

  iMat F, uF;
  sim::manifold_boundary_inner(T, F, uF);
  
  EXPECT_EQ(F.rows(), 6 * 2);
  EXPECT_EQ(uF.rows(), 4);
}

TEST(manifold_boundary_inner, edge_num)
{
  dMat V;
  iMat F;
  igl::readOBJ(SHARED_PATH "/3holes.obj", V, F);

  iMat E, uE;
  sim::manifold_boundary_inner(F, E, uE);

  EXPECT_EQ(E.rows(), 0);
  EXPECT_EQ( V.rows() - uE.rows() + F.rows(), 2 - 2 * 3);
}

TEST(manifold_boundary_inner, face_bound_attr)
{

  iMat F(2, 3);
  F << 0, 1, 2,
       3, 2, 1;
  iVec FA(2);
  FA << 0, 1;
  
  iMat E;
  iVec EA;
  sim::manifold_boundary_inner(F, FA, E, EA);

  EXPECT_EQ(E.rows(), 4);
  EXPECT_EQ(EA.sum(), 2);
}

TEST(manifold_boundary_inner, tet_bound_attr)
{
  /* tet cube
  *   4 - - - 7
  *   |\      |\
  *   | 5 - - - 6
  *   0 | - - 3 |
  *    \|      \|
  *     1 - - - 2
  * */
  iMat T(5, 4);
  T << 0, 1, 4, 3,
       2, 3, 6, 1,
       4, 5, 6, 1,
       4, 7, 6, 3,
       1, 3, 6, 4;

  for (int it = 0; it < REPEAT_N; it++)
  {
    dVec TA = dVec::Random(5);
    iMat F;
    dVec FA;
    sim::manifold_boundary_inner(T, TA, F, FA);

    EXPECT_EQ(F.rows(), 6 * 2);

    double sum = 3 * TA.head(4).sum();
    EXPECT_DOUBLE_EQ(FA.sum(), sum);
  }
}

TEST(manifold_boundary_inner, face_inner_attr)
{

  iMat F(2, 3);
  F << 0, 1, 2,
    3, 2, 1;
  iVec FA(2);
  FA << 0, 1;

  iMat E;
  iVec EA;
  iMat uE;
  iVec uEA1, uEA2;
  sim::manifold_boundary_inner(F, FA, E, EA, uE, uEA1, uEA2);

  EXPECT_EQ(E.rows(), 4);
  EXPECT_EQ(EA.sum(), 2);

  int uea = (uEA1 + uEA2).array().sum();
  EXPECT_EQ(uea, 1);
}

TEST(manifold_boundary_inner, tet_inner_attr)
{
  /* tet cube
  *   4 - - - 7
  *   |\      |\
  *   | 5 - - - 6
  *   0 | - - 3 |
  *    \|      \|
  *     1 - - - 2
  * */
  iMat T(5, 4);
  T << 0, 1, 4, 3,
    2, 3, 6, 1,
    4, 5, 6, 1,
    4, 7, 6, 3,
    1, 3, 6, 4;

  dVec TA(5);
  TA << 0, 1, 2, 3, 4;
  iMat F;
  dVec FA;
  iMat uF;
  dVec uFA1, uFA2;
  sim::manifold_boundary_inner(T, TA, F, FA, uF, uFA1, uFA2);

  EXPECT_EQ(F.rows(), 6 * 2);

  double sum = 3 * TA.head(4).sum();
  EXPECT_DOUBLE_EQ(FA.sum(), sum);

  dVec uFA = uFA1 + uFA2;
  dVec uFA_sort;
  iVec unuseI;
  igl::sort(uFA, 1, true, uFA_sort, unuseI);

  dVec uFA_gt(4);
  uFA_gt << 4, 5, 6, 7;

  EXPECT_EQ((uFA_gt - uFA_sort).norm(), 0);

}