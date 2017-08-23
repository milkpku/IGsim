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

#ifndef IGSIM_NEOHOOKEAN_MODEL_DPIOLA_H
#define IGSIM_NEOHOOKEAN_MODEL_DPIOLA_H
#include "igsim_inline.h"

#include <Eigen/Dense>

#include <vector>

namespace sim
{
  /*  Construct delta Piola tensor for given deformation gradient F, delta
   *  deformation gradient dF and parameter Mu, Lam
   *
   *  Neohookean model:
   *    dP(F, dF) = \mu * dF + [\mu - \lambda * log J] * F^{-T} dF^T F^{-T}
   *               + \lambda * tr(F^{-1} dF) F^{-T}
   *
   *  Inputs:
   *    F     3 by 3 matrix of deformation gradient
   *    dF    3 by 3 matrix of delta deformation gradient
   *    Mu    scalar of neohookean parameter \mu
   *    Lam   scalar of neohookean parameter \lambda
   *
   *  Outputs:
   *    dP    3 by 3 matrix of delta Piola tensor caused by dF
   */
  template <
    typename DerivedF, typename DeriveddF,
    typename ScalarMu, typename ScalarLam,
    typename DeriveddP>
  IGSIM_INLINE void neohookean_model_dPiola(
    const Eigen::PlainObjectBase<DerivedF>& F,
    const Eigen::PlainObjectBase<DeriveddF>& dF,
    const ScalarMu& Mu,
    const ScalarLam& Lam,
    Eigen::PlainObjectBase<DeriveddP>& dP);
  /*  Inputs:
   *    dF  std::vector of 3 by 3 matrix of delta deformation gradient
   *  
   *  Outputs:
   *    dP  std::vector of 3 by 3 matrix of delta Piola tensor caused by dF
   */
  template <
    typename DerivedF, typename DeriveddF_T, typename DeriveddF_A,
    typename ScalarMu, typename ScalarLam,
    typename DeriveddP_T, typename DeriveddP_A>
  IGSIM_INLINE void neohookean_model_dPiola(
    const Eigen::PlainObjectBase<DerivedF>& F,
    const std::vector<DeriveddF_T, DeriveddF_A>& dF,
    const ScalarMu& Mu,
    const ScalarLam& Lam,
    std::vector<DeriveddP_T, DeriveddP_A>& dP);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "neohookean_model_dPiola.cpp"
#endif


#endif //!IGSIM_NEOHOOKEAN_MODEL_DPIOLA_H
