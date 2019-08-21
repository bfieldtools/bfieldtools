import numpy as np
from scipy.sparse import csr_matrix, spdiags
from .utils import dual_areas


def laplacian_matrix(mesh):
    """
    Sparse Laplace(-Beltrami) operator

    Parameters:

        mesh: Trimesh Mesh object

    Returns:

        Cotangent weights: w_ij = - 0.5* (cot(alpha) + cot(beta))

    """
    N = mesh.vertices.shape[0]
    R = mesh.vertices[mesh.faces]  # Nt x 3 (corners) x 3 (xyz)

    # Edges opposite to the vertex
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2) # Nt x 3 (edges) x 3 (x,y,z)

    ii = []
    jj = []
    cot = []
    # Loop over edges in triangles
    for i in range(3):
        i1 = (i+1) % 3
        i2 = (i+2) % 3
        ii.append(mesh.faces[:, i1])
        jj.append(mesh.faces[:, i2])
        cot.append(-0.5*(edges[:, i1, :]*edges[:, i2, :]).sum(axis=-1)/(2*mesh.area_faces))

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

def mass_matrix(mesh, da=None):
    '''
    Computes mass matrix of mesh.

    Parameters:

        mesh: Trimesh Mesh object

    Returns:
        Mesh mass matrix (Nvertices, Nvertices)
    '''
    if da is None:
        da = dual_areas(mesh.faces, mesh.area_faces)

    A = spdiags(da, 0, mesh.vertices.shape[0], mesh.vertices.shape[0]).tocsr()

    return A

def gradient_matrix(mesh, rotated=False):
    """
    Calculate a (rotated) gradient matrix for hat basis functions
    (stream functions) in the triangular mesh described by

    Parameters:

        mesh: Trimesh mesh object
        rotated: Boolean argument describing whether the gradient should be rotated 90 degrees

    Returns:
        Gx ,Gy, Gx (Ntris, Nverts) matrices for calculating the components
        of gradient at triangles
    """
    R = mesh.vertices[mesh.faces]  # Nt x 3 (corners) x 3 (xyz)
    # Calculate edge vectors for each triangle
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)


    tri_data = edges/(2*mesh.area_faces[:, None, None])
    if not rotated:
        # Rotate 90 degrees CW to get the original gradient
        tri_data = np.cross(mesh.face_normals[:, None, :], tri_data, axis=-1)

    tri_data = tri_data.reshape(-1, 3).T
    ii = np.array([[i]*3 for i in range(len(mesh.faces))]).ravel()  # [0,0,0,1,1,1,...]
    jj = mesh.faces.ravel()  # [t[0,0], t[0,1], t[0,2], t[1,0], t[1,1], t[1,2], ...]
    Gx = csr_matrix((tri_data[0], (ii, jj)),
                    shape=(mesh.faces.shape[0], mesh.vertices.shape[0]), dtype=float)
    Gy = csr_matrix((tri_data[1], (ii, jj)),
                    shape=(mesh.faces.shape[0], mesh.vertices.shape[0]), dtype=float)
    Gz = csr_matrix((tri_data[2], (ii, jj)),
                    shape=(mesh.faces.shape[0], mesh.vertices.shape[0]), dtype=float)
    return Gx, Gy, Gz

def gradient(vals, mesh, rotated=False):
    """
    Applies mesh (rotated) gradient matrix operator on vector that is
    defined in the vertex locations of the mesh

    Parameters:
        vals: Nv x 1 array of data to compute the gradient of
        mesh: Trimesh object describing the triangular mesh
        rotated: Boolean argument describing whether the gradient should be rotated 90 degrees

    returns:
        gradient (3, Ntris)
    """
    Gx, Gy, Gz = gradient_matrix(mesh, rotated)
    return np.array([Gx @ vals, Gy @ vals, Gz @ vals])
