'''
Contains functions for calculating the coupling of surface current density in a
triangle mesh to magnetic field.

'''

import numpy as np
from numba import jit
import time

from .utils import get_quad_points
from .laplacian_mesh import gradient_matrix
from .integrals import triangle_potential_dipole_linear, triangle_potential_uniform


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


def create_basis(mesh, centre=np.array([0, 0, 0])):
    '''
    Calculate "hat" basis functions for each vertex in a given mesh, see Michael Poole's thesis.

    Parameters
    ----------
    mesh: Trimesh mesh object

    Returns
    -------
    basis: dict
        dict containing basis functions for each vertex
    '''
    n_verts = len(mesh.vertices)

    vmi = []
    A = []
    B = []
    C = []
    for m in range(n_verts):

        vmi.append([])

        A.append([])
        B.append([])
        C.append([])

        vert_links, vert_idx = np.where(mesh.faces == m)

        for tri_idx in vert_links:

            #indices within triangle
            m_vert_idx = np.where(mesh.faces[tri_idx] == m)[0][0]
            a_vert_idx = (m_vert_idx + 1) % 3
            b_vert_idx = (m_vert_idx + 2) % 3


            #A - coordinates of the second node of the i_th face of the m_th node
            Ami = mesh.vertices[mesh.faces[tri_idx]][a_vert_idx]
            #B - coordinates of the third node of the i_th face of the m_th node
            Bmi = mesh.vertices[mesh.faces[tri_idx]][b_vert_idx]
            #C - coordinates of the first node of the i_th face of the m_th node
            Cmi = mesh.vertices[m]

            A[m].append(Ami)
            B[m].append(Bmi)
            C[m].append(Cmi)


            vectorOA = Ami - centre

            dotprod = np.dot(mesh.face_normals[tri_idx], vectorOA)


            if dotprod > 0:
                vmi[m].append((Ami - Bmi)/(2 * mesh.area_faces[tri_idx]))
            else:
                vmi[m].append((Bmi - Ami)/(2 * mesh.area_faces[tri_idx]))

    return dict(v=vmi, A=A, B=B, C=C)


def compute_C_loops(mesh, r, basis=None, vert_links=None):
    '''
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.
    See eq. 5.13 in Michael Poole's thesis.

    Parameters
    -----------
    mesh: Trimesh mesh object describing mesh
    r: (Np, 3) array
        Field evaluation points
    basis: dict
        basis functions used in computation can be given as parameter
    vert_links: list of lists
        Describes the neighborhood for each vertex

    Returns
    -------
    C: (Np, Nvertices, 3) array
        Coupling matrix for surface current in the mesh to the evaluation points

    '''
    mu0 = 4 * np.pi * 1e-7
    coef = mu0 / (4 * np.pi)


    print('Computing C matrix, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    if vert_links is None:
        vert_links = get_vert_links(mesh.vertices, mesh.faces)

    if basis is None:
        basis = create_basis(mesh)


    w_quad, r_quad = get_quad_points(mesh.vertices, mesh.faces, method='centroid')
#    n_quad_points = len(w_quad)

    n_target_points = len(r)

    n_verts = len(mesh.vertices)

    C = np.zeros((n_target_points, n_verts, 3))

    #Convert nested list structures to numpy arrays, numba can't handle nested lists
    vert_links_arr, n_links = make_2D_array(vert_links)
    bval_arr, n_links = make_3D_array(basis['v'])

    C = _compute_C_loops(mesh.vertices,
                         vert_links_arr.astype(int),
                         n_links,
                         r,
                         mesh.area_faces,
                         r_quad,
                         w_quad,
                         bval_arr)

    duration = time.time() - start

    print('took %.2f seconds.'%duration)

    return coef * C


def make_2D_array(lis):
    """
    Function to get 2D array from a list of lists
    """
    n = len(lis)
    lengths = np.array([len(x) for x in lis])
    max_len = max(lengths)
    arr = np.zeros((n, max_len))

    for i in range(n):
        arr[i, :lengths[i]] = lis[i]
    return arr, lengths

def make_3D_array(lis):
    """
    Function to get 3D [x, y, 3] array from a list of lists of 3x1 vectors
    """
    n = len(lis)
    lengths = np.array([len(x) for x in lis])
    max_len = max(lengths)
    arr = np.zeros((n, max_len, 3))

    for i in range(n):
        arr[i, :lengths[i]] = lis[i]
    return arr, lengths


@jit(nopython=True, parallel=True, fastmath=True, nogil=True)
def _compute_C_loops(verts, vert_links, n_links, r, tri_areas, r_quad, w_quad, basis_value):
    '''
    C matrix computation backend, uses numba for speed and parallelization.
    '''

    n_target_points = len(r)
    n_verts = len(verts)
    n_quad_points = len(w_quad)

    C_part = np.zeros((n_target_points, n_verts, 3))

    #Initialize variables
    element = np.array([0., 0., 0.])
    denom = 0.

    #For each vertex
    for n in range(n_verts):

#        print('vertex: %d'%n)

        #For each target point
        for k in range(n_target_points):

            #For each triangle the vertex is used for
            for i in range(n_links[n]):

                element = np.array([0., 0., 0.])

                #For each quadrature point of that triangle
                for l in range(n_quad_points):
                    denom = np.linalg.norm(r[k] - r_quad[vert_links[n][i]][l])**3
#                    element += w_quad[l] * mycross(-(r[k] - r_quad[vert_links[n][i]][l]) / denom, basis_value[n][i]).flatten()

                    #Faster to do component-wise than using vector, numba doesn't support numpy cross product
                    element[0] += w_quad[l] * ((-basis_value[n][i][2]*(r[k, 1] - r_quad[vert_links[n][i]][l, 1]) + basis_value[n][i][1] * (r[k, 2] - r_quad[vert_links[n][i]][l, 2]))) / denom
                    element[1] += w_quad[l] * ((-basis_value[n][i][0]*(r[k, 2] - r_quad[vert_links[n][i]][l, 2]) + basis_value[n][i][2] * (r[k, 0] - r_quad[vert_links[n][i]][l, 0]))) / denom
                    element[2] += w_quad[l] * ((-basis_value[n][i][1]*(r[k, 0] - r_quad[vert_links[n][i]][l, 0]) + basis_value[n][i][0] * (r[k, 1] - r_quad[vert_links[n][i]][l, 1]))) / denom

                #Area integral
                C_part[k, n] += element * tri_areas[vert_links[n][i]]

    return C_part

def compute_C(mesh, r, Nchunks=None):
    '''
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object describing mesh
    r: target points (Np, 3)

    Returns
    -------
    C: (Np, 3, Nvertices) array
        Coupling matrix for surface current in the mesh to the evaluation points)

    '''
    mu0 = 4 * np.pi * 1e-7
    coef = mu0 / (4 * np.pi)

    print('Computing C matrix, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    w_quad, r_quad = get_quad_points(mesh.vertices, mesh.faces, method='centroid')

    # Rotated gradients (currents)
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)

    # Initialize C-matrix
    n_target_points = len(r)
    n_verts = len(mesh.vertices)
    C = np.zeros((n_target_points, n_verts, 3))

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

    for n in range(Nchunks):
        # Diffence vectors (Neval, Ntri, Nquad, 3)
        RR = r_quad[None, :, :, :] - r[n::Nchunks, None, None, :]

        # RR/norm(RR)**3 "Gradient of Green's function"
        g = - RR/((np.linalg.norm(RR, axis=-1)**3)[:, :, :, None])

        # Sum over quad points and multiply by triangle area
        g = (g*w_quad[:, None]).sum(axis=-2)
        g *= mesh.area_faces[:, None]

        # Cross product RR/norm(RR)
        C[n::Nchunks, :, 0] = g[:, :, 2] @ Gy - g[:, :, 1] @ Gz
        C[n::Nchunks, :, 1] = g[:, :, 0] @ Gz - g[:, :, 2] @ Gx
        C[n::Nchunks, :, 2] = g[:, :, 1] @ Gx - g[:, :, 0] @ Gy


    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    return coef * np.moveaxis(C, 2, 1)

def compute_C_analytic(mesh, r, Nchunks=None):
    '''
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object describing the mesh
    r: target points (Np, 3)

    Returns
    -------
    C: (Np, 3, Nvertices) array
        Coupling matrix for surface current in the mesh to the evaluation points)

    '''
    from .integrals import omega, gamma0
    coef = 1e-7

    print('Computing C matrix, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

    ta = mesh.area_faces
    tn = mesh.face_normals

    # Nfaces, 3, 3
    rfaces = mesh.vertices[mesh.faces]
    # Neval, Nfaces, Nedges
    coeffs = np.zeros((r.shape[0:1] + rfaces.shape[:-1]))
    # Nfaces, Nedges, 3
    edges = np.roll(rfaces, 1, -2) - np.roll(rfaces, 2, -2)
    grad = np.cross(tn[:, None, :], edges, axis=-1)/(2*ta[:, None, None])
    solid_angle = np.zeros((r.shape[0], rfaces.shape[0]))
    # Calculate potentials and related coefficients
    for n in range(Nchunks):
        RRchunk = r[n::Nchunks, None, None, :] - rfaces[None, :, :, :]
        # Neval, Nfaces, Nedges
        result = -np.sum(np.sum(gamma0(RRchunk)[..., None]*edges,
                               axis=-2)[...,None, :]*edges, axis=-1)
        result *= 1/(2*ta[..., :, None])
        solid_angle[n::Nchunks] = omega(RRchunk)
        coeffs[n::Nchunks] = result

#    # Accumulate the elements
    C = np.zeros((r.shape[0], mesh.vertices.shape[0], 3))
    for ind_f, f in enumerate(mesh.faces):
        C[:, f, :] += coeffs[:, ind_f, :, None]*tn[ind_f]
        C[:, f, :] += -solid_angle[:, ind_f:ind_f+1, None]*grad[ind_f]


    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    return coef * np.moveaxis(C, 2, 1)


def compute_U(mesh, r, Nchunks=None):
    """ Compute scalar potential matrix from linear stream functions
        using analytic integral
    """

    coeff = 1e-7  # mu_0/(4*pi)

    print('Computing U matrix, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    # Source and eval locations
    R1 = mesh.vertices[mesh.faces]
    R2 = r

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

    R2chunks = np.array_split(R2, Nchunks, axis=0)
    i0=0
    Uf = np.zeros((R2.shape[0], mesh.faces.shape[0], 3))
    for R2chunk in R2chunks:
        RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
        i1 = i0+RRchunk.shape[0]
        Pi = triangle_potential_dipole_linear(RRchunk, mesh.face_normals,
                                             mesh.area_faces)
        Uf[i0:i1] = Pi
        i0=i1

#     Accumulate the elements
    Uv = np.zeros((R2.shape[0], mesh.vertices.shape[0]))
    for ind_f, f in enumerate(mesh.faces):
        Uv[:, f] += Uf[:, ind_f]


    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    return Uv*coeff


def compute_A(mesh, r, Nchunks=None):
    """ Compute vector potential matrices (one for each coordinate)
        from linear stream functions using analytic integral
    """

    coeff = 1e-7  # mu_0/(4*pi)

    # Source and eval locations
    R1 = mesh.vertices[mesh.faces]
    R2 = r

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

    R2chunks = np.array_split(R2, Nchunks, axis=0)
    i0=0
    Af = np.zeros((R2.shape[0], mesh.faces.shape[0]))
    print('Computing potential matrix')

    for R2chunk in R2chunks:
        RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
        i1 = i0+RRchunk.shape[0]
        Pi = triangle_potential_uniform(RRchunk, mesh.face_normals, False)
        Af[i0:i1] = Pi
#        print((100*i1)//R2.shape[0], '% computed')
        i0=i1

    # Rotated gradients (currents)
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
    # Accumulate the elements
    Av = np.array([Af @ Gx, Af @ Gy, Af @ Gz])

    return Av*coeff
