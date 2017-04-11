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

#include "../include/IGsim/neohookean_model.h"
#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>

#define _USE_MATH_DEFINES
#include <math.h>

#define REPEAT_N 20
#define NUM_EQ(A, B) EXPECT_LT(abs(A-B), 1e-10)

TEST(neohookean_model, rigid_rotation)
{
  for(int i = 0; i < REPEAT_N; i++)
  {
    Eigen::Matrix3d F;
    F = Eigen::AngleAxisd(rand() * M_PI, Eigen::Vector3d::UnitX())
      * Eigen::AngleAxisd(rand() * M_PI, Eigen::Vector3d::UnitY())
      * Eigen::AngleAxisd(rand() * M_PI, Eigen::Vector3d::UnitZ());
    double mu = 1;
    double lambda = 1;

    double energy;
    Eigen::MatrixXd P, dPmu, dPlam;

    sim::neohookean_model(F, mu, lambda, energy);
	NUM_EQ(energy, 0.0);

    sim::neohookean_model(F, mu, lambda, energy, P);
	NUM_EQ(energy, 0.0);
	NUM_EQ(P.sum(), 0.0);
	NUM_EQ(P.rows(), 3);
	NUM_EQ(P.cols(), 3);

    sim::neohookean_model(F, mu, lambda, energy, P, dPmu, dPlam);
	NUM_EQ(energy, 0.0);
    NUM_EQ(P.sum(), 0.0);
	NUM_EQ(P.rows(), 3);
	NUM_EQ(P.cols(), 3);
	NUM_EQ(dPmu.sum(), 0.0);
	NUM_EQ(dPmu.rows(), 3);
	NUM_EQ(dPmu.cols(), 3);
	NUM_EQ(dPlam.sum(), 0.0);
	NUM_EQ(dPlam.rows(), 3);
	NUM_EQ(dPlam.cols(), 3);
  }
}

TEST(neohookean_model, increasing_energy)
{
  /* for tet */
  Eigen::Matrix3d F;
}

TEST(neohookean_model, increasing_force)
{
}

TEST(neohookean_model, increasing_material)
{
}
