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

#ifndef IGSIM_NEOHOOKEAN_MODEL_H
#define IGSIM_NEOHOOKEAN_MODEL_H
#include "igsim_inline.h"

#include <Eigen/Dense>

#include <vector>

namespace sim
{
  /*  Construct the energy scalar, Piola tensor, for given deformation gradient 
   *  F and material parameter Mu, Lam
   *
   *  Neohookean model: 
   *    energy(F) = \mu/2 * (I_1 - log I_3 - 3) + \lambda/8 * log^2 I_3 
   *    P(F) = \mu * (F - F^{-T}) + \lambda/2 * log I_3 * F^{-T}
   *    
   *    where 
   *    I_1 = tr(F^T * F), I_3 = det(F^T * F) = det^2(F) = J^2
   *
   *  Inputs:
   *    F     3 by 3 matrix of deformation gradient
   *    Mu    scalar of neohookean parameter \mu
   *    Lam   scalar of neohookean parameter \lambda 
   *
   *  Outputs:
   *    energy  scalar, deformation energy densisty
   */
  template < 
    typename DerivedF, typename ScalarMu, typename ScalarLam, 
    typename ScalarE>
  IGSIM_INLINE void neohookean_model(
    const Eigen::PlainObjectBase<DerivedF>& F,
    const ScalarMu& Mu,
    const ScalarLam& Lam,
    ScalarE& energy);
  /*  Outputs:
   *    P   3 by 3 matrix of Piola tensor    
   */
  template < 
    typename DerivedF, typename ScalarMu, typename ScalarLam, 
    typename ScalarE, typename DerivedP>
  IGSIM_INLINE void neohookean_model(
    const Eigen::PlainObjectBase<DerivedF>& F,
    const ScalarMu& Mu,
    const ScalarLam& Lam,
    ScalarE& energy,
    Eigen::PlainObjectBase<DerivedP>& P);
  /*  Outputs:
   *    dPmu  3 by 3 matrix of dPiola tensor by dmu
   *    dPlam 3 by 3 matrix of dPiola tensor by dlam
   */
  template < 
    typename DerivedF, typename ScalarMu, typename ScalarLam, 
    typename ScalarE, typename DerivedP, 
    typename DerivedPmu, typename DerivedPlam>
  IGSIM_INLINE void neohookean_model(
    const Eigen::PlainObjectBase<DerivedF>& F,
    const ScalarMu& Mu,
    const ScalarLam& Lam,
    ScalarE& energy,
    Eigen::PlainObjectBase<DerivedP>& P,
    Eigen::PlainObjectBase<DerivedPmu>& dPmu,
    Eigen::PlainObjectBase<DerivedPlam>& dPlam);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "neohookean_model.cpp"
#endif

#endif // !IGSIM_NEOHOOKEAN_MODEL_H
