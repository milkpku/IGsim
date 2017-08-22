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

#include "elastic_neohookean.h"
#include "neohookean_model.h"
#include "neohookean_model_dPiola.h"

#include "average_onto_faces_mat.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

template <
  typename DerivedV, typename DerivedT,
  typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
  typename DerivedMu, typename DerivedLam,
  typename DerivedE>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedT>& T, 
  const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
  const Eigen::PlainObjectBase<DerivedW>& W,
  const Eigen::PlainObjectBase<DerivedMu>& _Mu,
  const Eigen::PlainObjectBase<DerivedLam>& _Lam,
  DerivedE&  energy)
{
  using namespace std;
  assert(V.cols() == 3 && "Vertices dim must be 3");
  assert(T.cols() == 4 && "Tetra dim must be 4");
  assert(Bm.size() == T.rows() && "Bm size missmatch T rows");
  assert(W.size() == T.rows() && "W size missmatch T rows");
  assert(_Mu.size() == _Lam.size() && "Mu size missmatch Lam size");
  
  DerivedMu Mu;
  DerivedLam Lam;
  if (_Mu.size() == V.rows())
  {
    /* average parameter \mu and \lambda from vertices to tetrahedrons */
    Eigen::SparseMatrix<DerivedE> Proj;
    sim::average_onto_faces_mat(V, T, Proj);

    Mu = Proj * _Mu;
    Lam = Proj * _Lam;
  }
  else if (_Mu.size() == T.rows())
  {
    Mu = _Mu;
    Lam = _Lam;
  }
  else
  {
    assert(false && "Mu size missmatch V.rows() and T.rows()");
    return;
  }
  typedef Eigen::Matrix3d Mat3;

  energy = 0;

  /* loop for each tet */
  for (int i = 0; i < T.rows(); i++)
  {
    /* compute deformation pos Ds = [v0-v3, v1-v3, v2-v3] */
    Mat3 Ds_t;
    Ds_t << V.row(T(i, 0)), V.row(T(i, 1)), V.row(T(i, 2));
    Ds_t.rowwise() -= V.row(T(i, 3));

    /* compute deformation gradient F */
    Mat3 F;
    F = Ds_t.transpose() * Bm[i];

    /* call neohookean model to get energy */
    DerivedE tmp_e;
    sim::neohookean_model(F, Mu(i), Lam(i), tmp_e);

    /* accumulate energy */
    energy += W(i) * tmp_e;
  if (energy != energy) return;
  }
 
}

template <
  typename DerivedV, typename DerivedT, 
  typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
  typename DerivedMu, typename DerivedLam,
  typename DerivedE, typename DerivedF>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedT>& T, 
  const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
  const Eigen::PlainObjectBase<DerivedW>& W,
  const Eigen::PlainObjectBase<DerivedMu>& _Mu,
  const Eigen::PlainObjectBase<DerivedLam>& _Lam,
  DerivedE&  energy,
  Eigen::PlainObjectBase<DerivedF>& f)
{
  using namespace std;
  assert(V.cols() == 3 && "Vertices dim must be 3");
  assert(T.cols() == 4 && "Tetra dim must be 4");
  assert(Bm.size() == T.rows() && "Bm size missmatch T rows");
  assert(W.size() == T.rows() && "W size missmatch T rows");
  assert(Mu.size() == Lam.size() && "Mu size missmatch Lam size");
  
  DerivedMu Mu;
  DerivedLam Lam;
  if (_Mu.size() == V.rows())
  {
    /* average parameter \mu and \lambda from vertices to tetrahedrons */
  Eigen::SparseMatrix<DerivedE> Proj;
  sim::average_onto_faces_mat(V, T, Proj);

  Mu = Proj * _Mu;
  Lam = Proj * _Lam;
  }
  else if (_Mu.size() == T.rows())
  {
  Mu = _Mu;
  Lam = _Lam;
  }
  else
  {
  assert(false && "Mu size missmatch V.rows() and T.rows()");
  return;
  }
  typedef Eigen::Matrix3d Mat3;

  energy = 0;

  f = DerivedF::Zero(V.rows(), V.cols());

  /* loop for of each tet */
  for (int i = 0; i < T.rows(); i++)
  {
    /* compute deformation pos Ds = [v0-v3, v1-v3, v2-v3] */
    Mat3 Ds_t;
    Ds_t << V.row(T(i, 0)), V.row(T(i, 1)), V.row(T(i, 2));
    Ds_t.rowwise() -= V.row(T(i, 3));

    /* compute deformation gradient F */
    Mat3 F;
    F = Ds_t.transpose() * Bm[i];

    /* call neohookean model to get energy */
    DerivedE tmp_e;
    Mat3 P;
    sim::neohookean_model(F, Mu(i), Lam(i), tmp_e, P);
    
    /* accumulate energy */
    energy += W(i) * tmp_e;
  if (energy != energy) return;

    /* accumulate forces */
    Mat3 _H = - W(i) * P * Bm[i].transpose();
    Mat3 H = _H.transpose();  /* change to be row wise force */
    f.row(T(i, 0)) += H.row(0);
    f.row(T(i, 1)) += H.row(1);
    f.row(T(i, 2)) += H.row(2);
    f.row(T(i, 3)) -= H.colwise().sum();
  }

}

template <
  typename DerivedV, typename DerivedT, 
  typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
  typename DerivedMu, typename DerivedLam,
  typename DerivedE, typename DerivedF, typename ScalarK>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedT>& T, 
  const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
  const Eigen::PlainObjectBase<DerivedW>& W,
  const Eigen::PlainObjectBase<DerivedMu>& _Mu,
  const Eigen::PlainObjectBase<DerivedLam>& _Lam,
  DerivedE&  energy,
  Eigen::PlainObjectBase<DerivedF>& f,
  Eigen::SparseMatrix<ScalarK>& K)
{
  using namespace std;
  assert(V.cols() == 3 && "Vertices dim must be 3");
  assert(T.cols() == 4 && "Tetra dim must be 4");
  assert(Bm.size() == T.rows() && "Bm size missmatch T rows");
  assert(W.size() == T.rows() && "W size missmatch T rows");
  assert(_Mu.size() == _Lam.size() && "Mu size missmatch Lam size");
  
  DerivedMu Mu;
  DerivedLam Lam;
  if (_Mu.size() == V.rows())
  {
    /* average parameter \mu and \lambda from vertices to tetrahedrons */
    Eigen::SparseMatrix<DerivedE> Proj;
    sim::average_onto_faces_mat(V, T, Proj);

    Mu = Proj * _Mu;
    Lam = Proj * _Lam;
  }
  else if (_Mu.size() == T.rows())
  {
  Mu = _Mu;
  Lam = _Lam;
  }
  else
  {
  assert(false && "Mu size missmatch V.rows() and T.rows()");
  return;
  }
  typedef Eigen::Matrix3d Mat3;

  energy = 0;

  f = DerivedF::Zero(V.rows(), V.cols());

  K.setZero();
  K.resize(V.size(), V.size());
  typedef Eigen::Triplet<ScalarK> Tk;
  std::vector<Tk> K_coeff;
  K_coeff.clear();
  /* 12 parameters of a tet influence each other */
  K_coeff.reserve(144 * T.rows());

  /* loop for of each tet */
  for (int i = 0; i < T.rows(); i++)
  {
    /* compute deformation pos Ds = [v0-v3, v1-v3, v2-v3] */
    Mat3 Ds_t;
    Ds_t << V.row(T(i, 0)), V.row(T(i, 1)), V.row(T(i, 2));
    Ds_t.rowwise() -= V.row(T(i, 3));

    /* compute deformation gradient F */
    Mat3 F;
    F = Ds_t.transpose() * Bm[i];

    /* call neohookean model to get energy */
    DerivedE tmp_e;
    Mat3 P;
    sim::neohookean_model(F, Mu(i), Lam(i), tmp_e, P);
    
    /* accumulate energy */
    energy += W(i) * tmp_e;
  if (energy != energy) return;

    /* accumulate forces */
    Mat3 _H = - W(i) * P * Bm[i].transpose();
    Mat3 H = _H.transpose();  /* change to be row wise force */
    f.row(T(i, 0)) += H.row(0);
    f.row(T(i, 1)) += H.row(1);
    f.row(T(i, 2)) += H.row(2);
    f.row(T(i, 3)) -= H.colwise().sum();

    /* accumulate stiffness matrix */
    std::vector<Mat3, Eigen::aligned_allocator<Mat3>> dF, dP;
    dF.clear();
    dF.reserve(12);
    /* Conpute the influence of k th coordinate of j th vertex in tet i on
     * total forces by pipline dx -> dF -> dP -> dH -> df
     *  j is the index of vertex in tet, 
     *  k is the index of coordinate x, y, z 
     */
    /* compute dx -> dF, j = 0, 1, 2 */
    for(int j = 0; j < 3; j++)
      for(int k = 0; k < 3; k++)
      {
        Mat3 dDs = Eigen::Matrix3d::Zero();
        dDs(k, j) = 1;

        Mat3 dF_tmp = dDs * Bm[i];
        dF.push_back(dF_tmp);
      }
    /* compute dx -> dF j = 3, namely the last vertex in tet */
    for(int k = 0; k < 3; k++)
    { 
      Mat3 dDs = Eigen::Matrix3d::Zero();
      dDs.row(k) << -1, -1, -1;

      Mat3 dF_tmp = dDs * Bm[i];
      dF.push_back(dF_tmp);
    }

    /* batch compute dF -> dP */
    sim::neohookean_model_dPiola(F, dF, Mu(i), Lam(i), dP); 
    
    /* compute dP -> dH -> df
     * address of k th coordinate of jth vertex of tet is 
     * k * V.rows() + T(i, j) according to K's definition 
     */
    int count = 0;
    int num_V = V.rows();
    /* compute dP -> dH -> df, j= 0, 1, 2, 3, the index of source vertex 
     * (j, k) means the k th part of j th vertex of tet T.row(i)*/
    for(int j = 0; j < 4; j++)
      for(int k = 0; k < 3; k++)
      {
        Mat3 _dH = - W(i) * dP[count++] * Bm[i].transpose();
        Mat3 dH = _dH.transpose(); /* change to be row wise force */
        int src_addr = k * num_V + T(i, j);
        
        /* compute dH -> df, l = 0, 1, 2, the index of target vertex */
        for(int l = 0; l < 3; l++)
          for(int m = 0; m < 3; m++)
          {
            int targ_addr = m * num_V + T(i, l);
            K_coeff.push_back(Tk(targ_addr, src_addr, dH(l, m)));
          }
        /* compute dH -> df, l = 3 */
        auto df_3 = - dH.colwise().sum();
        for(int m = 0; m < 3; m++)
        {
          int targ_addr = m * num_V + T(i, 3);
          K_coeff.push_back(Tk(targ_addr, src_addr, df_3(m)));
        }
      }
  }
  K.setFromTriplets(K_coeff.begin(), K_coeff.end());
  K.makeCompressed();
}

template <
  typename DerivedV, typename DerivedT,
  typename DerivedBm_T, typename DerivedBm_A, typename DerivedW,
  typename DerivedMu, typename DerivedLam,
    typename DerivedFMu, typename DerivedFLam,
  typename ScalarMu, typename ScalarLam>
  IGSIM_INLINE void sim::elastic_neohookean(
    const Eigen::PlainObjectBase<DerivedV>& V,
    const Eigen::PlainObjectBase<DerivedT>& T,
    const std::vector<DerivedBm_T, DerivedBm_A>& Bm,
    const Eigen::PlainObjectBase<DerivedW>& W,
    const Eigen::PlainObjectBase<DerivedMu>& _Mu,
    const Eigen::PlainObjectBase<DerivedLam>& _Lam,
    Eigen::PlainObjectBase<DerivedFMu>& fmu,
    Eigen::PlainObjectBase<DerivedFLam>& flam,
    Eigen::SparseMatrix<ScalarMu>& Kmu,
    Eigen::SparseMatrix<ScalarLam>& Klam)
{
  using namespace std;
  assert(V.cols() == 3 && "Vertices dim must be 3");
  assert(T.cols() == 4 && "Tetra dim must be 4");
  assert(Bm.size() == T.rows() && "Bm size missmatch T rows");
  assert(W.size() == T.rows() && "W size missmatch T rows");
  assert(_Mu.size() == _Lam.size() && "Mu size missmatch Lam size");
  
  DerivedMu Mu;
  DerivedLam Lam;
  if (_Mu.size() == V.rows())
  {
    /* average parameter \mu and \lambda from vertices to tetrahedrons */
    Eigen::SparseMatrix<double> Proj;
    sim::average_onto_faces_mat(V, T, Proj);

    Mu = Proj * _Mu;
    Lam = Proj * _Lam;
  }
  else if (_Mu.size() == T.rows())
  {
    Mu = _Mu;
    Lam = _Lam;
  }
  else
  {
    assert(false && "Mu size missmatch V.rows() and T.rows()");
    return;
  }
  typedef Eigen::Matrix3d Mat3;

  /* initialize fmu and flam */
  fmu.resize(Mu.size());
  flam.resize(Lam.size());

  /* initialize Kmu and Klam */
  typedef Eigen::Triplet<ScalarMu> Tmu;
  Kmu.setZero();
  Kmu.resize(V.size(), Mu.size());
  std::vector<Tmu> Kmu_coeff;
  Kmu_coeff.clear();
  Kmu_coeff.reserve(12 * T.rows());

  typedef Eigen::Triplet<ScalarLam> Tlam;
  Klam.setZero();
  Klam.resize(V.size(), Lam.size());
  std::vector<Tlam> Klam_coeff;
  Klam_coeff.clear();
  Klam_coeff.reserve(12 * T.rows());

  /* loop for of each tet */
  for (int i = 0; i < T.rows(); i++)
  {
    /* compute deformation pos Ds = [v0-v3, v1-v3, v2-v3] */
    Mat3 Ds_t;
    Ds_t << V.row(T(i, 0)), V.row(T(i, 1)), V.row(T(i, 2));
    Ds_t.rowwise() -= V.row(T(i, 3));

    /* compute deformation gradient F */
    Mat3 F;
    F = Ds_t.transpose() * Bm[i];

    /* call neohookean model to get dPmu, dPlam */
    ScalarMu demu;
    ScalarLam delam;
    Mat3 dPmu, dPlam;
    sim::neohookean_model(F, Mu(i), Lam(i), demu, delam, dPmu, dPlam);
    int num_V = V.rows();

    /* fill in fmu and flam */
    fmu(i) = demu;
    flam(i) = delam;

    ////////////////////// fill in Kmu /////////////////////
    Mat3 _dHmu = - W(i) * dPmu * Bm[i].transpose();
    Mat3 dHmu = _dHmu.transpose();
    /* compute dH -> df, j = 0, 1, 2 */
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
        int targ_addr = k * num_V + T(i, j);
        Kmu_coeff.push_back(Tmu(targ_addr, i, dHmu(j, k)));
      }
    /* compute dH -> df, j = 3 */
    auto dfmu_3 = - dHmu.colwise().sum();
    for (int k = 0; k < 3; k++)
    {
      int targ_addr = k * num_V + T(i, 3);
      Kmu_coeff.push_back(Tmu(targ_addr, i, dfmu_3(k)));
    }
    
    ///////////////////// fill in Klam /////////////////////
    Mat3 _dHlam = - W(i) * dPlam * Bm[i].transpose();
    Mat3 dHlam = _dHlam.transpose();
    /* compute dH -> df, j = 0, 1, 2 */
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
      {
        int targ_addr = k * num_V + T(i, j);
        Klam_coeff.push_back(Tlam(targ_addr, i, dHlam(j, k)));
      }
    /* compute dH -> df, j = 3 */
    auto dflam_3 = - dHlam.colwise().sum();
    for (int k = 0; k < 3; k++)
    {
      int targ_addr = k * num_V + T(i, 3);
      Klam_coeff.push_back(Tlam(targ_addr, i, dflam_3(k)));
    }
  }

  Kmu.setFromTriplets(Kmu_coeff.begin(), Kmu_coeff.end());
  Klam.setFromTriplets(Klam_coeff.begin(), Klam_coeff.end());
}

template <
  typename DerivedV, typename DerivedT, 
  typename DerivedBm_T, typename DerivedBm_A, typename DerivedW>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedT>& T,
  std::vector<DerivedBm_T, DerivedBm_A>& Bm,
  Eigen::PlainObjectBase<DerivedW>& W)
{
  assert(V.cols() == 3 && "Vertices dim must be 3");
  assert(T.cols() == 4 && "Tetra dim must be 4");
  typedef Eigen::Matrix3d Mat3;
  
  Bm.clear();
  Bm.reserve(T.rows());
  W = DerivedW::Zero(T.rows());

  for(int i = 0; i < T.rows(); i++)
  {
    Mat3 Ds_t;
    Ds_t << V.row(T(i, 0)), V.row(T(i, 1)), V.row(T(i, 2));
    Ds_t.rowwise() -= V.row(T(i, 3));

    W(i) = Ds_t.determinant();
    DerivedBm_T Bm_tmp = Ds_t.inverse().transpose();
    
    Bm.push_back(Bm_tmp);
  }

  assert( W.minCoeff() > 0 && "negtive volume in W");
}

template <
  typename DerivedV, typename DerivedVinit, typename DerivedT,
  typename DerivedMu, typename DerivedLam,
  typename DerivedE>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedVinit>& Vinit,
  const Eigen::PlainObjectBase<DerivedT>& T, 
  const Eigen::PlainObjectBase<DerivedMu>& Mu,
  const Eigen::PlainObjectBase<DerivedLam>& Lam,
  DerivedE&  energy)
{
  typedef Eigen::Matrix3d Mat3;
  std::vector<Mat3, Eigen::aligned_allocator<Mat3>> Bm;
  Eigen::VectorXd W;

  /* precomputing Bm and W */
  elastic_neohookean(Vinit, T, Bm, W);

  elastic_neohookean(V, T, Bm, W, Mu, Lam, energy);
}

template <
  typename DerivedV, typename DerivedVinit, typename DerivedT,
  typename DerivedMu, typename DerivedLam,
  typename DerivedE, typename DerivedF>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedVinit>& Vinit,
  const Eigen::PlainObjectBase<DerivedT>& T, 
  const Eigen::PlainObjectBase<DerivedMu>& Mu,
  const Eigen::PlainObjectBase<DerivedLam>& Lam,
  DerivedE&  energy,
  Eigen::PlainObjectBase<DerivedF>& f)
{
  typedef Eigen::Matrix3d Mat3;
  std::vector<Mat3, Eigen::aligned_allocator<Mat3>> Bm;
  Eigen::VectorXd W;

  /* precomputing Bm and W */
  elastic_neohookean(Vinit, T, Bm, W);

  elastic_neohookean(V, T, Bm, W, Mu, Lam, energy, f);
}

template <
  typename DerivedV, typename DerivedVinit, typename DerivedT,
  typename DerivedMu, typename DerivedLam,
  typename DerivedE, typename DerivedF, typename ScalarK>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedVinit>& Vinit,
  const Eigen::PlainObjectBase<DerivedT>& T, 
  const Eigen::PlainObjectBase<DerivedMu>& Mu,
  const Eigen::PlainObjectBase<DerivedLam>& Lam,
  DerivedE&  energy,
  Eigen::PlainObjectBase<DerivedF>& f,
  Eigen::SparseMatrix<ScalarK>& K)
{
  typedef Eigen::Matrix3d Mat3;
  std::vector<Mat3, Eigen::aligned_allocator<Mat3>> Bm;
  Eigen::VectorXd W;

  /* precomputing Bm and W */
  elastic_neohookean(Vinit, T, Bm, W);

  elastic_neohookean(V, T, Bm, W, Mu, Lam, energy, f, K);
}

template <
  typename DerivedV, typename DerivedVinit, typename DerivedT,
  typename DerivedMu, typename DerivedLam,
  typename DerivedFMu, typename DerivedFLam,
  typename ScalarMu, typename ScalarLam>
IGSIM_INLINE void sim::elastic_neohookean(
  const Eigen::PlainObjectBase<DerivedV>& V,
  const Eigen::PlainObjectBase<DerivedVinit>& Vinit,
  const Eigen::PlainObjectBase<DerivedT>& T, 
  const Eigen::PlainObjectBase<DerivedMu>& Mu,
  const Eigen::PlainObjectBase<DerivedLam>& Lam,
  Eigen::PlainObjectBase<DerivedFMu>& fmu,
  Eigen::PlainObjectBase<DerivedFLam>& flam,
  Eigen::SparseMatrix<ScalarMu>& Kmu,
  Eigen::SparseMatrix<ScalarLam>& Klam)
{
  typedef Eigen::Matrix3d Mat3;
  std::vector<Mat3, Eigen::aligned_allocator<Mat3>> Bm;
  Eigen::VectorXd W;

  /* precomputing Bm and W */
  elastic_neohookean(Vinit, T, Bm, W);

  elastic_neohookean(V, T, Bm, W, Mu, Lam, fmu, flam, Kmu, Klam);
}

#ifdef IGSIM_STATIC_LIBRARY
#endif
