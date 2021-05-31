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

#ifndef IGSIM_ELASTIC_NEOHOOKEAN_H
#define IGSIM_ELASTIC_NEOHOOKEAN_H
#include "igsim_inline.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace sim
{
  /*  Precompute Bm and W for given tetmesh (Vinit, T)  for further computing
   *
   *  Inputs:
   *    V   #V by 3 matrix of vertices coordinate
   *    T   #T by 4 matrix of indices of tetrahedron corners into Vinit 
   *
   *  Outputs:
   *    Bm  #T std::vector of 3x3 matrix, inverse of [v0-v3, v1-v3, v2-v3]; 
   *    W   #T vector of volume of tetrahedron
   */
  template <
    typename DerivedV, typename DerivedT, 
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    Eigen::PlainObjectBase<DerivedW>& W);

  /*  Construct the enregy scalar, force vector,
   *  stiffness matrix for a given tetmesh (V, Bm, T)
   *
   *  Inputs:
   *    V   #V by 3 matrix of vertex coordinates
   *    T   #T by 4 matrix of indices of tetrahedron corners into V
   *    Bm  std::vector by #T of 3x3 matrix, inverse of [v0-v3, v1-v3, v2-v3]
   *    W   #T vector of weighting/volume of tetrahedron
   *    Mu  {#T|#V} vector of neohookean parameter \mu on tet/vertex
   *    Lam {#T|#V} vector of neohookean parameter \lambda on tet/vertex
   *
   *  Outputs:
   *    energy  scalar, total energy of neohookean elastic system
   */
  template <
    typename DerivedV, typename DerivedT,
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam,
    typename DerivedE>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T, 
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& Mu,
    const Eigen::PlainObjectBase<DerivedLam>& Lam,
    DerivedE&  energy);
  /*  Outputs:
   *    f   #V by 3 matrix of gradient of energy on vertices
   */
  template <
    typename DerivedV, typename DerivedT, 
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam,
    typename DerivedE, typename DerivedG>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T, 
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& Mu,
    const Eigen::PlainObjectBase<DerivedLam>& Lam,
    DerivedE&  energy,
    Eigen::PlainObjectBase<DerivedG>& grad);
  /*  Outputs:
   *    K   (3 * #V) by (3 * #V) sparse matrix of hessian matrix, 
   *      dgrad = H dx, where df and dx are (3 * #V) vectors of delta force and
   *      delta pos.
   *      dgrad << dgrad_x, dgrad_y, dgrad_z, 
   *      dx << dx_x, dx_y, dx_z, 
   *      where dgrad_* and dx_* are (#V) vectors of corresponding delta scalar.
   *
   *      This representation is convenient for ColMajor Eigen Matrix to reshape
   *      from #V by 3 matrix to (3 * #V) vector.
   */
  template <
    typename DerivedV, typename DerivedT, 
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam,
    typename DerivedE, typename DerivedG, typename ScalarH>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T, 
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& Mu,
    const Eigen::PlainObjectBase<DerivedLam>& Lam,
    DerivedE&  energy,
    Eigen::PlainObjectBase<DerivedG>& grad,
    Eigen::SparseMatrix<ScalarH>& H);
   /* Outputs:
    *   fmu  #Mu vector of gradient of energy on \mu
    *   flam #Lam vector of gradient of energy on \lambda
    *   Kmu  (3 * #V) by #Mu sparse matrix of gradiant of force on \mu
    *   Klam (3 * #V) by #Lam sparse matrix of gradiant of force on \lambda
    *      both Kmu and Klam correspond to df defined above.
    */
  template <
    typename DerivedV, typename DerivedT, 
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam,
    typename DerivedFMu, typename DerivedFLam,
    typename ScalarMu, typename ScalarLam>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T, 
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& Mu,
    const Eigen::PlainObjectBase<DerivedLam>& Lam,
    Eigen::PlainObjectBase<DerivedFMu>& fmu,
    Eigen::PlainObjectBase<DerivedFLam>& flam,
    Eigen::SparseMatrix<ScalarMu>& Kmu,
    Eigen::SparseMatrix<ScalarLam>& Klam);

  /* Only calculate properties for part of tets in T
   *
   * Inputs:
   *  ...
   *  start int, start row of T
   *  num   int, num of tets from start
   *  ...
   * */
  template <
    typename DerivedV, typename DerivedT,
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam,
    typename DerivedE>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T, 
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& Mu,
    const Eigen::PlainObjectBase<DerivedLam>& Lam,
    const int start,
    const int num,
    DerivedE&  energy);
  /*  Outputs:
   *    f   #V by 3 matrix of forces on vertices
   */
  template <
    typename DerivedV, typename DerivedT, 
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam,
    typename DerivedE, typename DerivedG>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& _Mu,
    const Eigen::PlainObjectBase<DerivedLam>& _Lam,
    const int start,
    const int num,
    DerivedE& energy,
    Eigen::PlainObjectBase<DerivedG>& grad);
  /*  Outputs:
   *    H   (3 * #V) by (3 * #V) sparse matrix of hessian matrix, 
   */
  template <
    typename DerivedV, typename DerivedT, 
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam,
    typename DerivedE, typename DerivedG, typename ScalarH>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& _Mu,
    const Eigen::PlainObjectBase<DerivedLam>& _Lam,
    const int start,
    const int num,
    DerivedE& energy,
    Eigen::PlainObjectBase<DerivedG>& grad,
    std::vector<Eigen::Triplet<ScalarH>>& hess_vec);

  /*  Inputs:
   *    basis   (3 * #V) by N dense matrix of basis
   *
   *  Outputs:
   *    hess    N by N dense matrix, hessian of energy based on N basis
   */
  template<
    typename DerivedV, typename DerivedT, 
    typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
    typename DerivedMu, typename DerivedLam, typename DerivedB,
    typename DerivedE, typename DerivedG, typename DerivedH>
  IGSIM_INLINE void elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& _Mu,
    const Eigen::PlainObjectBase<DerivedLam>& _Lam,
    const Eigen::PlainObjectBase<DerivedB>& basis,
    const int start,
    const int num,
    DerivedE& energy,
    Eigen::PlainObjectBase<DerivedG>& grad,
    Eigen::PlainObjectBase<DerivedH>& hess);
}

#ifndef IGSIM_STATIC_LIBRARY
#   include "elastic_neohookean.cpp"
#endif

#endif //! IGSIM_ELASTIC_NEOHOOKEAN_H
