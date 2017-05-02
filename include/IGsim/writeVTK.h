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

#ifndef IGSIM_WRITEVTK_H
#define IGSIM_WRITEVTK_H
#include "igsim_inline.h"

#include <Eigen/Dense>

#include <fstream>
#include <string>

namespace sim
{
  enum VTK_TYPE{VTK_SCALAR, VTK_COLOR, VTK_VECTOR, VTK_TEXTURE};
  /* write information of tetrahedron mesh to vtk file,
   * including vertex coordinate and tetrahedron connectivity
   *
   *  Inputs:
   *    filename string, output filename
   *    V #V by 3 matrix of vertices coordinate
   *    T #T by 4 matrix of indices of tetrahedron vertices into V
   */
  template <typename DerivedV, typename DerivedT>
  IGSIM_INLINE void writeVTK(
      const std::string& filename,
      const Eigen::PlainObjectBase<DerivedV>& V,
      const Eigen::PlainObjectBase<DerivedT>& T);
  /*
   *  also include attribute of vertices/tetrahedrons
   *
   *    info_type VTK_TYPE, specify information type
   *    V_info #V/#T by {1|3} matrix of vertices {scalar|vector} information
   */
  template <typename DerivedV, typename DerivedT, typename DerivedVi>
  IGSIM_INLINE void writeVTK(
      const std::string& filename,
      const Eigen::PlainObjectBase<DerivedV>& V,
      const Eigen::PlainObjectBase<DerivedT>& T,
      const VTK_TYPE& info_type,
      const Eigen::PlainObjectBase<DerivedVi>& V_info);
  /* write information of tetrahedron mesh to vtk filestream
   * 
   *  Inputs:
   *    fout  ofstream, output filestream
   *    V #V by 3 matrix of vertices coordinate
   *    T #T by 4 matrix of indices of tetrahedron vertices into V
   */
  template <typename DerivedV, typename DerivedT>
  IGSIM_INLINE void writeVTK(
      std::ofstream& fout,
      const Eigen::PlainObjectBase<DerivedV>& V,
      const Eigen::PlainObjectBase<DerivedT>& T);
  /*
   *  Inputs:
   *    name  string, name of output information
   *    V_info #V/#T by {1|3} matrix of vertices {scalar|vector} information
   */
  template <typename DerivedVi>
  IGSIM_INLINE void writeVTK(
      std::ofstream& fout,
      const int V_size,
      const int T_size,
      const std::string name,
      const VTK_TYPE& info_type,
      const Eigen::PlainObjectBase<DerivedVi>& V_info);
}

#ifndef IGSIM_STATIC_LIBRARAY
#   include "writeVTK.cpp"
#endif

#endif
