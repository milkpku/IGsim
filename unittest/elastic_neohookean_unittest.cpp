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

using namespace Eigen;
typedef Eigen::MatrixXd dMat;
typedef Eigen::VectorXd dVec;
typedef Eigen::MatrixXi iMat;
typedef Eigen::Matrix3d dMat3;
typedef Eigen::SparseMatrix<double> SpMat;

void init_tetmesh(dMat& V, iMat& T)
{
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
}

TEST(elastic_neohookean, energy)
{
  dMat V;
  iMat T;
  init_tetmesh(V, T);
  std::vector<dMat3> Bm;
  dVec W;
  sim::elastic_neohookean(V, T, Bm, W);

  for (int i = 0; i < REPEAT_N; i++)
  {
    dMat dV = 0.1 * dMat::Random(V.rows(), V.cols());
    dMat V_curr = V + dV;

    dVec mu = dVec::Random(T.rows());
    dVec lam = dVec::Random(T.rows());

    double e_curr;


    sim::elastic_neohookean(V_curr, T, Bm, W, mu, lam, e_curr);

    EXPECT_TRUE(std::isfinite(e_curr));
  }
}

TEST(elastic_neohookean, force)
{
  dMat V;
  iMat T;
  init_tetmesh(V, T);
  std::vector<dMat3> Bm;
  dVec W;
  sim::elastic_neohookean(V, T, Bm, W);

  for (int i = 0; i < REPEAT_N; i++)
  {
    dMat dV = 0.3 * dMat::Random(V.rows(), V.cols());
    dMat V_curr = V + dV;

    dMat dV2 = dMat::Random(V.rows(), V.cols());
    dV2 *= 1e-5 / dV2.norm();
    dMat V_targ = V_curr + dV2;

    dVec mu = dVec::Random(T.rows());
    dVec lam = dVec::Random(T.rows());

    double e_curr, e_targ;
    dMat f_curr;

    sim::elastic_neohookean(V_curr, T, Bm, W, mu, lam, e_curr, f_curr);
    sim::elastic_neohookean(V_targ, T, Bm, W, mu, lam, e_targ);

    double de = e_targ - e_curr;

    double de_exp = -(f_curr.cwiseProduct(dV2)).sum();

    double de_err = de - de_exp;

    EXPECT_LT(abs(de_err / de), 3e-4);
  }
}

TEST(elastic_neohookean, stiffness)
{
  dMat V;
  iMat T;
  init_tetmesh(V, T);
  std::vector<dMat3> Bm;
  dVec W;
  sim::elastic_neohookean(V, T, Bm, W);

  for (int i = 0; i < REPEAT_N; i++)
  {
    dMat dV = 0.1 * dMat::Random(V.rows(), V.cols());
    dMat V_curr = V + dV;

    dMat dV2 = dMat::Random(V.rows(), V.cols());
    dV2 *= 1e-5 / dV2.norm();
    dMat V_targ = V_curr + dV2;

    dVec mu = dVec::Random(T.rows());
    dVec lam = dVec::Random(T.rows());

    double e_curr, e_targ;
    dMat f_curr, f_targ;
    SpMat K_curr;
    sim::elastic_neohookean(V_curr, T, Bm, W, mu, lam, e_curr, f_curr, K_curr);
    sim::elastic_neohookean(V_targ, T, Bm, W, mu, lam, e_targ, f_targ);

    dMat df = f_targ - f_curr;
    df.resize(df.size(), 1);

    dV2.resize(dV2.size(), 1);
    dMat df_exp = K_curr * dV2;

    dVec df_err = df - df_exp;

    EXPECT_LT(df_err.norm() / df.norm(), 1e-5);
  }
}

TEST(elastic_neohookean, dparam)
{
  dMat V;
  iMat T;
  init_tetmesh(V, T);
  std::vector<dMat3> Bm;
  dVec W;
  sim::elastic_neohookean(V, T, Bm, W);

  for (int i = 0; i < REPEAT_N; i++)
  {
    dMat dV = 0.1 * dMat::Random(V.rows(), V.cols());
    dMat V_curr = V + dV;

    dVec mu = dVec::Random(T.rows());
    dVec lam = dVec::Random(T.rows());

    dVec dmu = dVec::Random(T.rows());
    dVec dlam = dVec::Random(T.rows());

    dVec mu_targ = mu + dmu;
    dVec lam_targ = lam + dlam;

    double e_curr, e_targ;
    dMat F_curr, F_targ;

    sim::elastic_neohookean(V_curr, T, Bm, W, mu, lam, e_curr, F_curr);

    sim::elastic_neohookean(V_curr, T, Bm, W, mu_targ, lam_targ, e_targ, F_targ);

    double de = e_targ - e_curr;
    dMat dF = F_targ - F_curr;
    dF.resize(dF.size(), 1);
    dVec dF_vec = dF.col(0);

    dVec fmu, flam;
    SpMat Kmu, Klam;
    sim::elastic_neohookean(V_curr, T, Bm, W, mu, lam, fmu, flam, Kmu, Klam);

    double de_est = fmu.dot(dmu) + flam.dot(dlam);
    dVec dF_est = Kmu * dmu + Klam * dlam;

    NUM_EQ(de_est, de);
    dVec dF_err = dF_vec - dF_est;
    NUM_EQ(dF_err.norm(), 0);
  }

}
