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

#include "unittest_defines.h"

TEST(average_onto_faces_mat, test) {

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  //TODO read obj file
  Eigen::SparseMatrix<double> A;

  sim::average_onto_faces_mat(V, F, A);

  for (int i = 0; i < REPEAT_N; ++i) {
    auto V_val = Eigen::VectorXd::Random(V.rows());
    auto F_weight = Eigen::VectorXd::Random(F.rows());

    double val1 = F_weight.dot(A * V_val);

    double val2 = 0;
    for (int j = 0; j < F.rows(); ++j) {
      auto f = F.row(i);
      double e = F_weight(j) * (V_val(f(0)) + V_val(f(1)) + V_val(f(2)));
      val2 += e / 3;
    }

    NUM_EQ(val1, val2);
  }

}
