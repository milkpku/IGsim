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

#include "../include/IGsim/elastic_neohookean.h"
#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "unittest_defines.h"

// Since elastic_neohookean E, force and K are tested by other means
// only test parameter here

TEST(elastic_neohookean, dparam)
{
  using namespace Eigen;
  MatrixXd V;
  MatrixXi T;
  
  // unit cube at origin
  V.resize(8, 3);
  T.resize(6, 4);
  V << 0, 0, 0,
       1, 0, 0,
       1, 1, 0,
       0, 1, 0,
       0, 0, 1,
       1, 0, 1,
       1, 1, 1,
       0, 1, 1;
  T << 1, 2, 5, 0,
       2, 6, 5, 0,
       5, 6, 4, 0,
       6, 7, 4, 0,
       6, 2, 7, 0,
       2, 3, 7, 0;

  for (int i = 0; i < REPEAT_N; i++)
  {
    MatrixXd dV = 0.1 * MatrixXd::Random(V.rows(), V.cols());
    MatrixXd V_curr = V + dV;

    VectorXd mu = VectorXd::Random(T.rows());
    VectorXd lam = VectorXd::Random(T.rows());

    VectorXd dmu = VectorXd::Random(T.rows());
    VectorXd dlam = VectorXd::Random(T.rows());

    VectorXd mu_targ = mu + dmu;
    VectorXd lam_targ = lam + dlam;

    double e_curr, e_targ;
    MatrixXd F_curr, F_targ;

    sim::elastic_neohookean(V_curr, V, T, mu, lam, e_curr, F_curr);

    sim::elastic_neohookean(V_curr, V, T, mu_targ, lam_targ, e_targ, F_targ);

    double de = e_targ - e_curr;
    MatrixXd dF = F_targ - F_curr;
    dF.resize(dF.size(), 1);
    VectorXd dF_vec = dF.col(0);

    VectorXd fmu, flam;
    SparseMatrix<double> Kmu, Klam;
    sim::elastic_neohookean(V_curr, V, T, mu, lam, fmu, flam, Kmu, Klam);

    double de_est = fmu.dot(dmu) + flam.dot(dlam);
    VectorXd dF_est = Kmu * dmu + Klam * dlam;

    NUM_EQ(de_est, de);
    VectorXd dF_err = dF_vec - dF_est;
    NUM_EQ(dF_err.norm(), 0);
  }

}
