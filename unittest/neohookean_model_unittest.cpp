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
#include "../include/IGsim/neohookean_model_dPiola.h"
#include "gtest/gtest.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>


#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

#include "unittest_defines.h"

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
  /* regular tet */
  Eigen::Matrix3d F;
  F << 1, 0, 0,
       0, 1, 0,
       0, 0, 1;

  double mu = 1;
  double lambda = 1;
  int N = 100;
  double step = 1e-1;

  double energy;
  sim::neohookean_model(F, mu, lambda, energy);

  /* test stretch energy */
  for(int i = 1; i < N; i++)
  {
    double next_energy;
    F(0,0) = 1 + i * step;
    sim::neohookean_model(F, mu, lambda, next_energy);
    EXPECT_GT(next_energy, energy);
    energy = next_energy;
  }

  F(0,0) = 1;
  sim::neohookean_model(F, mu, lambda, energy);
  /* test press energy */
  for( int i = 2; i < N; i++)
  {
    double next_energy;
    F(0,0) = 1.0 / i;
    sim::neohookean_model(F, mu, lambda, next_energy);
    EXPECT_GT(next_energy, energy);
    energy = next_energy;
  }

}

TEST(neohookean_model, Piola_energy_relation)
{
  double mu = 1;
  double lambda = 1;
  double scale = 3;
  int N = 100;
  double test_ratio = 1e-6;

  /* test stretch energy and force */
  for(int i = 1; i < N; i++)
  {
    /* random tet */
    Eigen::Matrix3d F = scale * Eigen::Matrix3d::Random();
    while(F.determinant() < 0)
    {
      F = scale * Eigen::Matrix3d::Random();
    }

    double energy;
    sim::neohookean_model(F, mu, lambda, energy);

    Eigen::MatrixXd P;
    sim::neohookean_model(F, mu, lambda, energy, P);
    for(int j = 0; j <3; j++)
      for(int k = 0; k < 3; k++)
      {
        double test_step = F(j, k) * test_ratio;
        Eigen::Matrix3d dF = Eigen::Matrix3d::Zero();
        dF(j,k) = test_step;

        double test_energy;
        F += dF;
        sim::neohookean_model(F, mu, lambda, test_energy);
        F -= dF;

        Eigen::MatrixXd dP;
        sim::neohookean_model_dPiola(F, dF, mu, lambda, dP);

        double error = test_energy - energy - P(j, k) * test_step
                       - .5 * dP(j,k) * test_step;
        double error_ratio = error / (test_step * energy);
        EXPECT_LT(abs(error_ratio), 10 * test_ratio);
      }
  }
}

TEST(neohookean_model, increasing_material)
{
  double scale = 3;
  int N = 100;
  double test_ratio = 1e-6;

  srand(time(0));

  /* test stretch energy and force */
  for(int i = 1; i < N; i++)
  {
    /* random tet */
    Eigen::Matrix3d F = scale * Eigen::Matrix3d::Random();
    while(F.determinant() < 0)
    {
      F = scale * Eigen::Matrix3d::Random();
    }

    /* random material */
    double mu = scale * rand() / RAND_MAX;
    double lambda= scale * rand() / RAND_MAX;
    double energy;
    sim::neohookean_model(F, mu, lambda, energy);

    Eigen::MatrixXd P, dPmu, dPlam;
    sim::neohookean_model(F, mu, lambda, energy, P, dPmu, dPlam);

    /* test dPmu */
    for(int j = 0; j <3; j++)
      for(int k = 0; k < 3; k++)
      {
        double test_step = mu * test_ratio;

        double test_energy;
        Eigen::MatrixXd test_P;
        sim::neohookean_model(F, mu + test_step, lambda, test_energy, test_P);

        double error = (test_P - P - test_step * dPmu).squaredNorm();
        double base = (test_P - P).squaredNorm();

        double error_ratio = error / base;
        EXPECT_LT(abs(error_ratio), 1e-10);
      }

    /* test dPlam */
    for(int j = 0; j <3; j++)
      for(int k = 0; k < 3; k++)
      {
        double test_step = lambda * test_ratio;

        double test_energy;
        Eigen::MatrixXd test_P;
        sim::neohookean_model(F, mu, lambda + test_step, test_energy, test_P);

        double error = (test_P - P - test_step * dPlam).squaredNorm();
        double base = (test_P - P).squaredNorm();

        double error_ratio = error / base;
        EXPECT_LT(abs(error_ratio), 1e-10);
      }

  }
}
