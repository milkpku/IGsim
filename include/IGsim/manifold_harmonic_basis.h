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

#ifndef IGSIM_MANIFOLD_HARMONIC_BASIS_H
#define IGSIM_MANIFOLD_HARMONIC_BASIS_H
#include "igsim_inline.h"

#include <Eigen/Dense>

namespace sim
{
  /*  Constrct the manifold harmonic bases for a manifold (V, F)
   *  Eigenvalues are orded by increasing abs() value
   *
   *  Inputs:
   *    V #V by dim matrix of vertex coordinates
   *    F #F by {3|4} matrix of indices of simplex vertex into V
   *    k int, num of bases needed
   *
   *  Outputs:
   *    S k complex vector of eigenvalues from small abs() to large
   *    U #V by k complex matrix of eigenvectors
   */
  template <typename DerivedV, typename DerivedF>
  IGSIM_INLINE void manifold_harmonic_basis(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    const int& k,
    Eigen::VectorXcd& S,
    Eigen::MatrixXcd& U);
  /*
   *  Inputs:
   *    sigma double, band center, from which the algorithm find closest 
   *          eigenvalues and corresponding eigenvectors
  */
  template <typename DerivedV, typename DerivedF>
  IGSIM_INLINE void manifold_harmonic_basis(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedF>& F,
    const int& k,
    const double& sigma,
    Eigen::VectorXcd& S,
    Eigen::MatrixXcd& U);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "manifold_harmonic_basis.cpp"
#endif

#endif //!IGSIM_MANIFOLD_HARMONIC_BASIS_H
