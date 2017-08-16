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

#include "../include/IGsim/average_onto_faces_mat.h"
#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <igl/readOBJ.h>

#include "unittest_defines.h"

TEST(average_onto_faces_mat, test) {
  using namespace Eigen;
  MatrixXd V;
  MatrixXi F;

  //read obj file
  igl::readOBJ(SHARED_PATH "/frog.obj", V, F);
  SparseMatrix<double> A;

  sim::average_onto_faces_mat(V, F, A);

  for (int i = 0; i < REPEAT_N; ++i) {
    VectorXd V_val = VectorXd::Random(V.rows());
    VectorXd F_avg = A * V_val;

    for (int j = 0; j < F.rows(); ++j) {
      auto f = F.row(j);
      double avg = 0;
      for (int k = 0; k < F.cols(); ++k) {
        avg += V_val(f(k));
      }
      avg /= F.cols();
      NUM_EQ(avg, F_avg(j));
    }
  }

}
