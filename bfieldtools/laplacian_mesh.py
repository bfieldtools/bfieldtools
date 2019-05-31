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

def gradient_matrix(verts, tris, tri_normals=None, tri_areas=None, rotated=False):
    """ Calculate a (rotated) gradient matrix for hat basis functions
        (stream functions) in the triangular mesh described by

        verts: Nv x 3 array of mesh vertices (coordinates)
        tris: Nt x 3 array of mesh triangles (indices to verts array)

        return:
            Gx ,Gy, Gx (Ntris, Nverts) matrices for calculating the components
            of gradient at triangles
    """
    R = verts[tris]  # Nt x 3 (corners) x 3 (xyz)
    # Calculate edge vectors for each triangle
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_normals, type(None)) or isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)

    tri_data = edges/(2*tri_areas[:, None, None])
    if not rotated:
        # Rotate 90 degrees CW to get the original gradient
        tri_data = np.cross(tri_normals[:,None,:], tri_data, axis=-1)

    tri_data = tri_data.reshape(-1, 3).T
    ii = np.array([[i]*3 for i in range(len(tris))]).ravel()  # [0,0,0,1,1,1,...]
    jj = tris.ravel()  # [t[0,0], t[0,1], t[0,2], t[1,0], t[1,1], t[1,2], ...]
    Gx = csr_matrix((tri_data[0], (ii, jj)),
                    shape=(tris.shape[0], verts.shape[0]), dtype=float)
    Gy = csr_matrix((tri_data[1], (ii, jj)),
                    shape=(tris.shape[0], verts.shape[0]), dtype=float)
    Gz = csr_matrix((tri_data[2], (ii, jj)),
                    shape=(tris.shape[0], verts.shape[0]), dtype=float)
    return Gx, Gy, Gz

def gradient(vals, verts, tris, tri_normals=None, tri_areas=None, rotated=False):
    """ Calculate a (rotated) gradient for function values in hat basis
        (stream functions) in the triangular mesh described by

        verts: Nv x 3 array of mesh vertices (coordinates)
        tris: Nt x 3 array of mesh triangles (indices to verts array)

        return:
            gradient (3, Ntris)
    """
    Gx, Gy, Gz = gradient_matrix(verts, tris, tri_normals, tri_areas, rotated)
    return np.array([Gx @ vals, Gy @ vals, Gz @ vals])



