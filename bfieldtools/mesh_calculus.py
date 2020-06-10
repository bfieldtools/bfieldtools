"""
Contains functions for computing vector calculus quantities on triangle surface meshes.
These include the gradient, rotated gradient, divergence, curl, Laplacian and mass matrices.
"""

__all__ = [
    "curl",
    "curl_matrix",
    "divergence",
    "divergence_matrix",
    "gradient",
    "gradient_matrix",
    "laplacian_matrix",
    "mass_matrix",
]

import numpy as np
from scipy.sparse import csr_matrix, spdiags


def laplacian_matrix(mesh, material_param=None):
    """
    Sparse Laplace-Beltrami operator.

    This function returns the so-called cotangent-Laplacian L, which can be 
    used to solve the Laplacian eigenvectors as L @ v_i = -k_i^2 M @ v_i, 
    where M is the mass matrix of the mesh (see mass_matrix).
    
    For meshes with boundaries, the matrix corresponds to zero-Neumann
    boundary conditions (the natural boundary condtion).
    
    When solving something with zero-Dirichlet condition
    only the inner vertices should be free parameters and
    the associated elements of the matrix should be used.
    
    For more discussion about the boundary conditions related
    to stream functions of divergence-free surface currents,
    see the associated publication about the physics (Makinen 2020).

    Parameters
    ----------
    mesh: Trimesh Mesh object
    material_param: array-like with length N_triangles
        material parameter for each triangle, used for resistance matrix

    Returns
    -------
    sparse csr_matrix (N_vertices, N_vertices)
        non-diagonal elements L_ij = - 0.5* (cot(alpha) + cot(beta))
        diagonal elements L_ii = - sum_{i!=j} L_ij 

    """
    if material_param is None:
        p = 1
    else:
        p = material_param
    N = mesh.vertices.shape[0]
    R = mesh.vertices[mesh.faces]  # Nt x 3 (corners) x 3 (xyz)

    # Edges opposite to the vertex
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)

    ii = []
    jj = []
    cot = []
    # Loop over edges in triangles
    for i in range(3):
        i1 = (i + 1) % 3
        i2 = (i + 2) % 3
        ii.append(mesh.faces[:, i1])
        jj.append(mesh.faces[:, i2])
        c = (
            -0.5
            * (edges[:, i1, :] * edges[:, i2, :]).sum(axis=-1)
            / (2 * mesh.area_faces)
        )
        # Append cot with cotangent terms multiplied by the material parameter
        cot.append(c * p)

    ii = np.ravel(ii)
    jj = np.ravel(jj)
    cot = np.ravel(cot)
    # Build the sparse matrix
    L = csr_matrix((cot, (ii, jj)), shape=(N, N), dtype=float)
    # Sum contribution from both triangles (alpha and beta angles)
    # neighbouring the edge
    L = L + L.T
    # Add the diagonal
    L = L - spdiags(L.sum(axis=0), 0, N, N)

    return L


def mass_matrix(mesh, lumped=False):
    """
    Computes mass matrix of a mesh.

    Parameters
    ----------
    mesh: Trimesh Mesh object
    lumped: Boolean
        If True, use lumped approximation of mass matrix. If False (default),
        compute exact matrix. See Reuter et al 2009, page 3 (DOI: 10.1016/j.cag.2009.03.005)

    Returns
    -------
    sparse csr_matrix of variable shape (see the description of laplacian_matrix)
        Mass matrix

    """

    if lumped:
        from .utils import dual_areas

        da = dual_areas(mesh.faces, mesh.area_faces)
        M = spdiags(da, 0, mesh.vertices.shape[0], mesh.vertices.shape[0]).tocsr()
    else:
        N = mesh.vertices.shape[0]

        ii = []
        jj = []
        area = []
        # Loop over edges in triangles
        for i in range(3):
            i1 = (i + 1) % 3
            i2 = (i + 2) % 3
            ii.append(mesh.faces[:, i1])
            jj.append(mesh.faces[:, i2])
            area.append(mesh.area_faces / 12)

        ii = np.ravel(ii)
        jj = np.ravel(jj)
        area = np.ravel(area)
        # Build sparse matrix
        M = csr_matrix((area, (ii, jj)), shape=(N, N), dtype=float)
        # Sum contribution from both triangles (alpha and beta angles)
        # neighbouring the edge
        M = M + M.T
        M = M + spdiags(M.sum(axis=0), 0, N, N)

    return M


def gradient_matrix(mesh, rotated=False):
    """
    Calculate a (rotated) gradient matrix for hat basis functions
    (stream functions) in a mesh.

    Parameters
    ----------
    mesh: Trimesh mesh object
    rotated: boolean
        If True, rotate gradient 90 degrees clockwise

    Returns
    -------
    3 sparse csr_matrices
    Gx ,Gy, Gx (Ntris, Nverts) matrices
        for calculating the components of gradient at triangles
    """
    R = mesh.vertices[mesh.faces]  # Nt x 3 (corners) x 3 (xyz)
    # Calculate edge vectors for each triangle
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)

    tri_data = edges / (2 * mesh.area_faces[:, None, None])
    if not rotated:
        # Rotate 90 degrees CW to get the original gradient
        tri_data = np.cross(mesh.face_normals[:, None, :], tri_data, axis=-1)

    tri_data = tri_data.reshape(-1, 3).T
    ii = np.array(
        [[i] * 3 for i in range(len(mesh.faces))]
    ).ravel()  # [0,0,0,1,1,1,...]
    jj = mesh.faces.ravel()  # [t[0,0], t[0,1], t[0,2], t[1,0], t[1,1], t[1,2], ...]
    Gx = csr_matrix(
        (tri_data[0], (ii, jj)),
        shape=(mesh.faces.shape[0], mesh.vertices.shape[0]),
        dtype=float,
    )
    Gy = csr_matrix(
        (tri_data[1], (ii, jj)),
        shape=(mesh.faces.shape[0], mesh.vertices.shape[0]),
        dtype=float,
    )
    Gz = csr_matrix(
        (tri_data[2], (ii, jj)),
        shape=(mesh.faces.shape[0], mesh.vertices.shape[0]),
        dtype=float,
    )
    return Gx, Gy, Gz


def gradient(vals, mesh, rotated=False):
    """
    Applies mesh (rotated) gradient matrix operator on vector that is
    defined in the vertex locations of the mesh.

    Parameters
    ----------
    vals: Nv x 1 array of scalar data to compute the gradient of
    mesh: Trimesh object describing the triangular mesh
    rotated: boolean
        If True, rotate gradient 90 degrees clockwise

    Returns
    -------
    ndarray (3, Ntris)
        surface gradient of vals on each triangle

    """
    Gx, Gy, Gz = gradient_matrix(mesh, rotated)
    return np.array([Gx @ vals, Gy @ vals, Gz @ vals])


def divergence_matrix(mesh):
    """ 
    Divergence of tangential vector field on mesh faces as a linear mapping.

    Parameters
    ----------
    mesh: Trimesh object

    Returns
    -------
    3 sparse csc_matrices
        Dx ,Dy, Dz (Nverts, Ntri) matrices

    """
    from .utils import tri_normals_and_areas

    Gx, Gy, Gz = gradient_matrix(mesh, rotated=False)
    n, a = tri_normals_and_areas(mesh.vertices, mesh.faces)
    A = spdiags(a, 0, a.shape[0], a.shape[0])
    return -Gx.T * A, -Gy.T * A, -Gz.T * A


def curl_matrix(mesh):
    """ 
    Adjoint curl of tangential vector field.

    Parameters
    ----------
    mesh: Trimesh object

    Returns
    -------
    3 sparse csc_matrices
        Cx ,Cy, Cz (Nverts, Ntri) matrices

    """
    from .utils import tri_normals_and_areas

    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
    n, a = tri_normals_and_areas(mesh.vertices, mesh.faces)
    A = spdiags(a, 0, a.shape[0], a.shape[0])
    return -Gx.T * A, -Gy.T * A, -Gz.T * A


def divergence(vecs, mesh):
    """ 
    Divergence mapping applied to tangential vector field 'vecs'

    Parameters
    ----------

    vecs: ndarray (3, Nfaces)
            vector field at mesh faces

    Returns
    -------
    ndarray (Nverts,)
        Divergence applied on the vector field

    """
    Dx, Dy, Dz = divergence_matrix(mesh)
    return Dx @ vecs[:, 0] + Dy @ vecs[:, 1] + Dz @ vecs[:, 2]


def curl(vecs, mesh):
    """ 
    Curl applied to tangential vector field

    Parameters
    ----------
    mesh: Trimesh object

    Returns
    -------
    ndarray (Nverts,)
        Curl applied on the vector field


    """
    Cx, Cy, Cz = curl_matrix(mesh)
    return Cx @ vecs[:, 0] + Cy @ vecs[:, 1] + Cz @ vecs[:, 2]
