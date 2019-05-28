import numpy as np
from scipy.sparse import csr_matrix, spdiags
from .utils import tri_normals_and_areas, dual_areas


def laplacian_matrix(verts, tris, tri_normals=None, tri_areas=None):
    """
    Sparse Laplace(-Beltrami) operator

    Parameters:

        verts: (n, 3) array (float)
        tris: (m, 3) array (int) - indices into the verts array


    Returns:

        Cotangent weights: w_ij = - 0.5* (cot(alpha) + cot(beta))

    """
    N = verts.shape[0]
    R = verts[tris]  # Nt x 3 (corners) x 3 (xyz)

    # Edges opposite to the vertex
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2) # Nt x 3 (edges) x 3 (x,y,z)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_normals, type(None)) or isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)
    ii = []
    jj = []
    cot = []
    # Loop over edges in triangles
    for i in range(3):
        i1 = (i+1) % 3
        i2 = (i+2) % 3
        ii.append(tris[:, i1])
        jj.append(tris[:, i2])
        cot.append(-0.5*(edges[:, i1, :]*edges[:, i2, :]).sum(axis=-1)/(2*tri_areas))

    ii = np.ravel(ii)
    jj = np.ravel(jj)
    cot = np.ravel(cot)
    # Build sparse matrix
    L = csr_matrix((cot, (ii, jj)), shape=(N, N), dtype=float)
    # Sum contribution from both triangles (alpha and beta angles)
    # neighbouring the edge
    L = L + L.T
    L = L - spdiags(L.sum(axis=0), 0, N, N)

    return L

def mass_matrix(verts, tris, tri_areas=None, da=None):
    '''
    Computes mass matrix of mesh.
    '''

    if da is None:
        if tri_areas is None:
            tri_normals, tri_areas = tri_normals_and_areas(verts, tris)

        da = dual_areas(tris, tri_areas)

    A =spdiags(da, 0, verts.shape[0], verts.shape[0]).tocsr()

    return A

