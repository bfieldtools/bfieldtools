import numpy as np
from numba import jit
import time

from .utils import get_quad_points


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

    Parameters:
        mesh: Trimesh mesh object
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


def compute_C(mesh, r, basis=None, vert_links=None):
    '''
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.
    See eq. 5.13 in Michael Poole's thesis.

    Parameters:
        mesh: Trimesh mesh object describing mesh
        r: target points (N, 3)
        basis: basis functions used in computation can be given as parameter (dict)
        vert_links: list of lists describing the neighborhood for each vertex
    '''
    mu0 = 4 * np.pi * 1e-7
    coef = mu0 / (4 * np.pi)


    print('Computing C matrix, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    if vert_links is None:
        vert_links = get_vert_links(mesh.vertices, mesh.faces)

    if basis is None:
        basis = create_basis(mesh)


    w_quad, r_quad = get_quad_points(mesh.vertices, mesh.faces, method='Centroid')
#    n_quad_points = len(w_quad)

    n_target_points = len(r)

    n_verts = len(mesh.vertices)

    C = np.zeros((n_target_points, n_verts, 3))

    #Convert nested list structures to numpy arrays, numba can't handle nested lists
    vert_links_arr, n_links = make_2D_array(vert_links)
    bval_arr, n_links = make_3D_array(basis['v'])

    C = _compute_C(mesh.vertices,
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
def _compute_C(verts, vert_links, n_links, r, tri_areas, r_quad, w_quad, basis_value):
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
