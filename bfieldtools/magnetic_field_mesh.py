import numpy as np
from numba import njit, jit, prange

from joblib import Parallel, delayed
import multiprocessing

from utils import tri_normals_and_areas, get_quad_points




def get_neighbour_vertices(vertices, edges):
    '''
    Computes and returns the neighbor vertex indices for all vertices
    '''
#    return [edge[edge[0] == index] for edge in edges if index in edge]
    vi = []
    for vert_idx in range(len(vertices)):
        vi.append([])
        for edge in edges:
            if edge[0] == vert_idx:
                vi[vert_idx].append(edge[1])
            elif edge[1] == vert_idx:
                vi[vert_idx].append(edge[0])

        vi[vert_idx] = np.unique(vi[vert_idx])

    return vi


def get_vert_links(verts, tris):
    '''
    Computes and returns the triangles that each vertex corresponds to.
    '''

    vert_links = []
    for m in range(len(verts)):
        vert_links.append(np.where(tris == m)[0])

    return vert_links


def create_basis(verts, tris, centre=np.array([0, 0, 0]), tri_areas=None, tri_normals=None):
    '''
    Calculate basis function for each triangle
    '''
    n_verts = len(verts)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_normals, type(None)) or isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)

#    vert_links = get_neighbour_vertices(verts, edges)

    vmi = []
    A = []
    B = []
    C = []
    for m in range(n_verts):

        vmi.append([])

        A.append([])
        B.append([])
        C.append([])

        vert_links, vert_idx = np.where(tris == m)

        for tri_idx in vert_links:

            #indices within triangle
            m_vert_idx = np.where(tris[tri_idx] == m)[0][0]
            a_vert_idx = (m_vert_idx + 1) % 3
            b_vert_idx = (m_vert_idx + 2) % 3


            #A - coordinates of the second node of the i_th face of the m_th node
            Ami = verts[tris[tri_idx]][a_vert_idx]
            #B - coordinates of the third node of the i_th face of the m_th node
            Bmi = verts[tris[tri_idx]][b_vert_idx]
            #C - coordinates of the first node of the i_th face of the m_th node
            Cmi = verts[m]

            A[m].append(Ami)
            B[m].append(Bmi)
            C[m].append(Cmi)


            vectorOA = Ami - centre

            dotprod = np.dot(tri_normals[tri_idx], vectorOA)


            if dotprod > 0:
                vmi[m].append((Ami - Bmi)/(2 * tri_areas[tri_idx]))
            else:
                vmi[m].append((Bmi - Ami)/(2 * tri_areas[tri_idx]))

    return dict(v=vmi, A=A, B=B, C=C)


def compute_L(verts, tris, basis=None,
              vert_links=None,
              tri_areas=None,
              mu0=4*np.pi*1e-7,
              parallel=True):
    '''
    Compute and return mutual inductance matrix. 
    NB! This is not properly tested, use the implementation in mutual_inductance_mesh
    '''

    #Compute vertex links if not provided
    if vert_links is None:
        vert_links = get_vert_links(verts, tris)


    #Compute triangle areas if not provided
    if tri_areas is None:
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)


    #Compute basis if not provided
    if basis is None:
        basis = create_basis(verts, tris, tri_areas=tri_areas)

    # Gauss Legendre integration
    w_quad, r_quad = get_quad_points(verts, tris, method='Centroid', index=1)
    n_quad_points = len(w_quad)

    coef = mu0 / (4*np.pi)

    n_verts = len(verts)

    L = np.zeros((n_verts, n_verts))

    # If specified, split L matrix computation into chunks done in parallel
    if parallel:
        #Use all available cores
        num_cores = multiprocessing.cpu_count()

        #Determine the chunking of C
        vert_ranges = []
        for i in range(num_cores):
            vert_ranges.append((int(n_verts/num_cores * i), int(n_verts/num_cores * (i+1))))

        #Compute in parallel
        L_parts = Parallel(n_jobs=num_cores,
                           max_nbytes=2e9,
                           verbose=51)(delayed(_L_part)(verts,
                                                        vert_range,
                                                        vert_links,
                                                        tri_areas,
                                                        r_quad,
                                                        w_quad,
                                                        basis['v'],
                                                        basis['A'],
                                                        basis['B'],
                                                        basis['C']
                                                        ) for vert_range in vert_ranges)

        print(L_parts[0].shape)
        #Assemble into whole matrix
        for i in range(len(L_parts)):
            L[vert_ranges[i][0]:vert_ranges[i][1], :] = L_parts[i]

    #Otherwise, compute L matrix in one piece, one thread
    else:
        L = _L_part(verts,
                   (0, len(verts)),
                   vert_links,
                   tri_areas,
                   r_quad,
                   w_quad,
                   basis['v'],
                   basis['A'],
                   basis['B'],
                   basis['C'])
#
#    #Fill in lower triangular matrix, L is symmetric
#    i_lower = np.tril_indices(L.shape[0], -1)
#    L[i_lower] = L.T[i_lower]

    return coef * L


def _L_part(verts, vert_range, vert_links, tri_areas, r_quad, w_quad, basis_value, basis_A, basis_B, basis_C):
    '''
    Computes part of inductance matrix, for parallelization
    '''

    n_verts = len(verts)
    n_quad_points = len(w_quad)

    L_part = np.zeros((vert_range[1] - vert_range[0], n_verts))
        # Iterate though upper triangle of matrix
    for m in range(vert_range[0], vert_range[1]):
        for n in range(m, n_verts):

            for i in range(len(vert_links[m])):

                currentTriangle_i = vert_links[m][i]

                for j in range(len(vert_links[n])):

                    currentTriangle_j = vert_links[n][j]

                    integral = 0

                    #If the current triangles are not the same,
                    #then we will not share any basis function. We can thus solve that easily
                    if currentTriangle_i != currentTriangle_j:
                        for k in range(n_quad_points):
                            for l in range(n_quad_points):
                                integral += w_quad[l]/np.linalg.norm(r_quad[currentTriangle_i, k] - r_quad[currentTriangle_j, l])

                            integral *= w_quad[k]

                        integral *= 2 * tri_areas[currentTriangle_j] *\
                                    2 * tri_areas[currentTriangle_i]

                    # If the 2 triangles are similar, we use an approximate
                    # calculation, according to page 72 equ 5.40 of Poole's
                    # thesis

                    else:
                        Ami = basis_A[m][i]
                        Bmi = basis_B[m][i]
                        Cmi = basis_C[m][i]

                        a = np.dot(Cmi - Ami, Cmi - Ami)
                        b = np.dot(Cmi - Ami, Cmi - Bmi)
                        c = np.dot(Cmi - Bmi, Cmi - Bmi)

                        sa = np.sqrt(a)
                        sc = np.sqrt(c)
                        ss = np.sqrt(a - 2*b + c)
                        sac = np.sqrt(a * c)

                        integral = (1 * (4 * tri_areas[currentTriangle_i]**2))                   \
                        * (1 / (6 * sa) * np.log(((a - b + sa * ss) * (b + sac)) /              \
                                               ((-b + sac) * (-a + b +sa*ss)))                  \
                        + 1 / (6 * sc) * np.log(((b + sac) * (-b + c + sc * ss)) /              \
                                                ((b - c + sc * ss) * (-b + sac)))               \
                        + 1 / (6 * ss) * np.log(((a - b + sa * ss) * (-b + c + sc * ss)) /      \
                                                ((b - c + sc * ss) * (-a + b + sa * ss) )))

                    L_part[m, n] += np.dot(basis_value[m][i], basis_value[n][j]) * integral

    return L_part

@jit
def compute_R(verts, tris, basis=None, vert_links=None, tri_areas=None, rho=1.68*1e-8, t=1e-4):
    '''
    Computes resitivity matrix for surface mesh made of material with resistivity rho and thickness t

    Equation come from the thesis of Michael Poole "Improved Equipment and
    %Techniques for Dynamic Shimming in High Field MRI
    page 66
    and G. N. Peeren, Stream function approach for determining optimal surface
    currents


    '''

    n_verts = len(verts)
    #Initialize sparse matrix
    R = np.full((n_verts, n_verts), fill_value=np.nan)



    surface_resistance = rho/t

    #Compute vertex links if needed
    if vert_links is None:
        vert_links = get_vert_links(verts, tris)

    #Compute basis functions if needded
    if basis is None:
        basis = create_basis(verts, tris)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)


    # Iterate though upper triangular matrix including diagonal
    for m in range(n_verts):
        for n in range(m, n_verts):
            bigSum = 0
            for i in range(len(vert_links[m])):

                currentTriangle_m = vert_links[m][i]


                for j in range(len(vert_links[n])):

                    currentTriangle_n = vert_links[n][j]


                    #if share a common node
                    if currentTriangle_m == currentTriangle_n:

                        vmi = basis['v'][m][i]

                        vnj = basis['v'][n][j]

                        bigSum += np.dot(vmi, vnj) * tri_areas[m]**2 #Square or not?

            R[m, n] = surface_resistance * bigSum


    #Fill in lower triangular matrix, R is symmetric
    i_lower = np.tril_indices(R.shape[0], -1)
    R[i_lower] = R.T[i_lower]

    return R


def compute_C(mesh, r, basis=None, vert_links=None, tri_areas=None, parallel=True):
    '''
    Given a mesh, computes the "cn matrix" or coupling matrix, see eq. 5.13 in Michael Poole's thesis.

    '''
    mu0 = 4 * np.pi * 1e-7
    coef = mu0 / (4 * np.pi)


    print('Computing Cn matrix...')

    if vert_links is None:
        vert_links = get_vert_links(mesh.vertices, mesh.faces)

    if tri_areas is None:
        tri_normals, tri_areas = tri_normals_and_areas(mesh.vertices, mesh.faces)

    if basis is None:
        basis = create_basis(mesh.vertices, mesh.faces)


    w_quad, r_quad = get_quad_points(mesh.vertices, mesh.faces, method='Centroid')
#    n_quad_points = len(w_quad)

    n_target_points = len(r)

    n_verts = len(mesh.vertices)

    C = np.zeros((n_target_points, n_verts, 3))

    # If specified, split C matrix computation into chunks done in parallel
    if parallel:

        #Use all available cores
        num_cores = multiprocessing.cpu_count()

        #Determine the chunking of C
        vert_ranges = []
        for i in range(num_cores):
            vert_ranges.append((int(n_verts/num_cores * i), int(n_verts/num_cores * (i+1))))

        #Compute in parallel
        Cn_parts = Parallel(n_jobs=num_cores,
                            max_nbytes=1e9,
                            verbose=51)(delayed(_Cn_part)(
                                       mesh.vertices[vert_range[0]:vert_range[1]],
                                       vert_links[vert_range[0]:vert_range[1]],
                                       r,
                                       tri_areas,
                                       r_quad,
                                       w_quad,
                                       basis['v'][vert_range[0]:vert_range[1]]
                                       ) for vert_range in vert_ranges)

        #Assemble into whole matrix
        for i in range(len(Cn_parts)):
            C[:, vert_ranges[i][0]:vert_ranges[i][1]] = Cn_parts[i]

    else:
        C = _Cn_part(mesh.vertices,
                      vert_links,
                      r,
                      tri_areas,
                      r_quad,
                      w_quad,
                      basis['v'])

    return coef * C


def _Cn_part(verts, vert_links, r, tri_areas, r_quad, w_quad, basis_value):
    '''
    Computes part of Cn matrix, used for parallelization.
    '''

    n_target_points = len(r)
    n_verts = len(verts)
    n_quad_points = len(w_quad)

    C_part = np.zeros((n_target_points, n_verts, 3))

    #For each vertex
    for n in range(n_verts):

#        print('vertex: %d'%n)

        #For each target point
        for k in range(n_target_points):

            #For each triangle the vertex is used for
            for i in range(len(vert_links[n])):

                element = 0.

                #For each quadrature point of that triangle
                for l in range(n_quad_points):
                    denom = np.linalg.norm(r[k] - r_quad[vert_links[n][i]][l])**3
                    element += w_quad[l] * np.cross(-(r[k] - r_quad[vert_links[n][i]][l]) / denom, basis_value[n][i]).flatten()

#                    C[k, n, 0] = w_quad[l] * ((-basis['v'][n][i][2]*(r[k, 1] - r_quad[vert_links[n][i]][l, 1]) + basis['v'][n][i][1] * (r[k, 2] - r_quad[vert_links[n][i]][l, 2]))) / denom
#                    C[k, n, 1] = w_quad[l] * ((-basis['v'][n][i][0]*(r[k, 2] - r_quad[vert_links[n][i]][l, 2]) + basis['v'][n][i][2] * (r[k, 0] - r_quad[vert_links[n][i]][l, 0]))) / denom
#                    C[k, n, 2] = w_quad[l] * ((-basis['v'][n][i][1]*(r[k, 0] - r_quad[vert_links[n][i]][l, 0]) + basis['v'][n][i][0] * (r[k, 1] - r_quad[vert_links[n][i]][l, 1]))) / denom

                #Area integral
                C_part[k, n] += element * tri_areas[vert_links[n][i]]

    return C_part
