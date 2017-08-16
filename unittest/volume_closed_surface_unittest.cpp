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

#include "../include/IGsim/volume_closed_surface.h"
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <igl/readOBJ.h>

#include "unittest_defines.h"

TEST(volume_closed_surface, tetrahedron) {

  using namespace Eigen;

  MatrixXd V;
  MatrixXi F;

  // tetrahedron, corner by plane (x + y + z = 1) and x-o-y, x-o-z, y-o-z plane 
  igl::readOBJ(SHARED_PATH "/tetrahedron.obj", V, F);

  // calculate volume, should be 1.0/6.0
  double vol;
  sim::volume_closed_surface(V, F, vol);
  NUM_EQ(vol, 1.0/6.0);
}

TEST(volume_closed_surface, cube) {

  using namespace Eigen;

  MatrixXd V;
  MatrixXi F;

  // unit cube
  igl::readOBJ(SHARED_PATH "/cube.obj", V, F);
  
  // calculate volume, should be 1.0
  double vol;
  sim::volume_closed_surface(V, F, vol);
  NUM_EQ(vol, 1.0);
}
