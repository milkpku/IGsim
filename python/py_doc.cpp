const char* __doc_sim_average_onto_faces_mat = R"sim_Qu8mg5v7(/*  Construct the average matrix from vertices to faces, which turn nodal
*  scalar to facial scalar
*
*  Inputs:
*    V #V by x matrix of vertex coordinates
*    F #F by x matrix, indices of simplex vertex into V
*
*  Outputs:
*    Proj  #F by #V sparse matrix, SF = Proj * SV, where SF and SV are scalar
*      represented on faces and vertices respectively
*/)sim_Qu8mg5v7";
const char* __doc_sim_boundary_facets = R"sim_Qu8mg5v7(/*  Construct boundary facet/line for tetrahedron/triangular meshes, and
*  transfer attribute attached on tetrahedron/triangular to facet/line
*
*  Inputs:
*    T #T by {3|4} matrix of ordered tetrahedron/triangualr mesh
*
*  Outputs:
*    F #F by {2|3} matrix of ordered face/line, boundary of given T
*/
/*
*  Inputs:
*    TA #T by n matrix of information for tetrahedron/triangular mesh
*
*  Outputs:
*    FA #F by n matrix of information for boundary face/line
*/)sim_Qu8mg5v7";
const char* __doc_sim_elastic_neohookean = R"sim_Qu8mg5v7(/*  Construct the enregy scalar, force vector,
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
/*  Outputs:
*    f   #V by 3 matrix of forces on vertices
*/
/*  Outputs:
*    K   (3 * #V) by (3 * #V) sparse matrix of force stiffness matrix,
*      df = K dx, where df and dx are (3 * #V) vectors of delta force and
*      delta pos.
*      df << df_x, df_y, df_z,
*      dx << dx_x, dx_y, dx_z,
*      where df_* and dx_* are (#V) vectors of corresponding delta scalar.
*
*      This representation is convenient for ColMajo Eigen Matrix to reshape
*      from #V by 3 matrix to (3 * #V) vector.
*/
/* Outputs:
*    Kmu  (3 * #V) by #Mu sparse matrix of gradiant of force on \mu
*    Klam (3 * #V) by #Lam sparse matrix of gradiant of force on \lambda
*      both Kmu and Klam correspond to df defined above.
*/
/*  Construct the enregy scalar, force vector,
*  stiffness matrix for a given tetmesh (V, Vinit, T)
*
*  Inputs:
*    V   #V by 3 matrix of vertex coordinates
*    Vinit #V by 3 matrix of vertex coordinates in rest pos
*    T   #T by 4 matrix of indices of tetrahedron corners into V
*    Mu  {#T|#V} vector of neohookean parameter \mu on tet/vertex
*    Lam {#T|#V} vector of neohookean parameter \lambda on tet/vertex
*
*  Outputs:
*    energy  scalar, total energy of neohookean elastic system
*/
/*  Outputs:
*    f   #V by 3 matrix of forces on vertices
*/
/*  Outputs:
*    K   (3 * #V) by (3 * #V) sparse matrix of force stiffness matrix,
*/
/* Outputs:
*    Kmu  (3 * #V) by #Mu sparse matrix of gradiant of force on \mu
*    Klam (3 * #V) by #Lam sparse matrix of gradiant of force on \lambda
*/
/*  Precompute Bm and W for given tetmesh (Vinit, T)
*
*  Inputs:
*    V   #V by 3 matrix of vertices coordinate
*    T   #T by 4 matrix of indices of tetrahedron corners into Vinit
*
*  Outputs:
*    Bm  #T std::vector of 3x3 matrix, inverse of [v0-v3, v1-v3, v2-v3];
*    W   #T vector of volume of tetrahedron
*/)sim_Qu8mg5v7";
const char* __doc_sim_igsim_inline = R"sim_Qu8mg5v7()sim_Qu8mg5v7";
const char* __doc_sim_manifold_harmonic_basis = R"sim_Qu8mg5v7(/*  Constrct the manifold harmonic bases for a manifold (V, F)
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
/*
*  Inputs:
*    sigma double, band center, from which the algorithm find closest
*          eigenvalues and corresponding eigenvectors
*/)sim_Qu8mg5v7";
const char* __doc_sim_manifold_volume = R"sim_Qu8mg5v7(/*  Calculate the manifold volume for each simplex (facet area for each
*  triangle in 2D manifold and volume for each tetrahedron in 3D manifold)
*
*  Inputs:
*    V #V by dim matrix of vertex coordinates
*    F #F by {3|4} matrix of indices of simplex vertex into V
*
*  Output:
*    Vol #F by 1 matrix of manifold volume
*/)sim_Qu8mg5v7";
const char* __doc_sim_neohookean_model = R"sim_Qu8mg5v7(/*  Construct the energy scalar, Piola tensor, for given deformation gradient
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
/*  Outputs:
*    P   3 by 3 matrix of Piola tensor
*/
/*  Outputs:
*    dPmu  3 by 3 matrix of dPiola tensor by dmu
*    dPlam 3 by 3 matrix of dPiola tensor by dlam
*/)sim_Qu8mg5v7";
const char* __doc_sim_neohookean_model_dPiola = R"sim_Qu8mg5v7(/*  Construct delta Piola tensor for given deformation gradient F, delta
*  deformation gradient dF and parameter Mu, Lam
*
*  Neohookean model:
*    dP(F, dF) = \mu * dF + [\mu - \lambda * log J] * F^{-T} dF F^{-T}
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
/*  Inputs:
*    dF  std::vector of 3 by 3 matrix of delta deformation gradient
*
*  Outputs:
*    dP  std::vector of 3 by 3 matrix of delta Piola tensor caused by dF
*/)sim_Qu8mg5v7";
const char* __doc_sim_tet_topology = R"sim_Qu8mg5v7(/*  Construct the topology of tetrahedrons
*
*  Inputs:
*    T   #T by 4 matrix of indices of vertices of tetrahedron corners,
*        ordered by [v1, v2, v3, v4], and the boudanry aim towards outside.
*
*  Outputs:
*    L   #T by 4 matrix of indices of neighbor tetrahedrons correspond to
*        vertices [v1, v2, v3, v4], or -1 when no neighbor tetrahedron
*/)sim_Qu8mg5v7";
const char* __doc_sim_tet_volume = R"sim_Qu8mg5v7(/*  Compute tetrahedron volume for given tetmesh, each tetrahedron is
*  [t0, t1, t2, t3] and its boundary towards outside is
*  [t0, t1, t2] - [t1, t2, t3] + [t2, t3, t0] - [t3, t0, t1]
*
*  Inputs:
*    V #V by 3 matrix of coordinate of vertices
*    T #T by 4 matrix of tetrahedron mesh
*
*  Output:
*    Vol #T vector of tetrahedron volume, with direction
*/)sim_Qu8mg5v7";
const char* __doc_sim_writeVTK = R"sim_Qu8mg5v7(/* write information of tetrahedron mesh to vtk file,
* including vertex coordinate and tetrahedron connectivity
*
*  Inputs:
*    filename string, output filename
*    V #V by 3 matrix of vertices coordinate
*    T #T by 4 matrix of indices of tetrahedron vertices into V
*/
/*
*  also include attribute of vertices/tetrahedrons
*
*    info_type VTK_TYPE, specify information type
*    V_info #V/#T by {1|3} matrix of vertices {scalar|vector} information
*/
/* write information of tetrahedron mesh to vtk filestream
*
*  Inputs:
*    fout  ofstream, output filestream
*    V #V by 3 matrix of vertices coordinate
*    T #T by 4 matrix of indices of tetrahedron vertices into V
*/
/*
*  Inputs:
*    name  string, name of output information
*    V_info #V/#T by {1|3} matrix of vertices {scalar|vector} information
*/)sim_Qu8mg5v7";
