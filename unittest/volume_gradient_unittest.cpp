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

#include "../include/IGsim/volume_gradient.h"
#include "../include/IGsim/volume_closed_surface.h"
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/readOBJ.h>
#include <stdlib.h>
#include <time.h>

#include "unittest_defines.h"

TEST(volume_gradient, gradient){
  using namespace std;
  using namespace Eigen;
  MatrixXd V;
  MatrixXi F;

  //read obj file
  igl::readOBJ(SHARED_PATH "/3holes.obj", V, F);

  // compute volume gradient
  double vol_start;
  sim::volume_closed_surface(V, F, vol_start);
  MatrixXd G;
  sim::volume_gradient(V, F, G);

  // randomly pick one node and test gradient 
  srand(time(NULL));
  for (int i = 0; i < REPEAT_N; ++i)
  {
    int id = rand() % V.rows();
    
    RowVectorXd dr = RowVectorXd::Random(V.cols());
    
    // compute volume change by volume gradient
    double dVol = dr.dot(G.row(id));

    // compute volume change by direct calculate volume
    V.row(id) += dr;
    double vol_end;
    sim::volume_closed_surface(V, F, vol_end);
    V.row(id) -= dr;

    NUM_EQ(vol_end - vol_start, dVol);
  }
}

TEST(volume_gradient, hessian){
  using namespace std;
  using namespace Eigen;
  MatrixXd V;
  MatrixXi F;

  //read obj file
  igl::readOBJ(SHARED_PATH "/3holes.obj", V, F);

  // compute volume hessian
  MatrixXd G;
  sim::volume_gradient(V, F, G);
  SparseMatrix<double> H;
  sim::volume_gradient(V, F, H);

  // randomly pick one node and test hessian 
  srand(time(NULL));
  for (int i = 0; i < REPEAT_N; ++i)
  {
    int id = rand() % V.rows();
    
    RowVectorXd dr = RowVectorXd::Random(V.cols());
    MatrixXd dV = MatrixXd::Zero(V.rows(), V.cols());
    dV.row(id) = dr;
    dV.resize(dV.size(), 1);
    
    // compute volume gradient change by volume hessian
    MatrixXd dG = H * dV;
    dG.resize(G.rows(), G.cols());

    // compute volume gradient change by direct calculate
    V.row(id) += dr;
    MatrixXd G_end;
    sim::volume_gradient(V, F, G_end);
    V.row(id) -= dr;

    MatrixXd diffG = G_end - G - dG;

    NUM_EQ(diffG.norm(), 0);
  }
}
