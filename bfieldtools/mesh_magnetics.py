"""
Contains functions for calculating the coupling of surface current density in a
triangle mesh to magnetic field as well as scalar and vector potentials.

"""

__all__ = [
    "magnetic_field_coupling",
    "magnetic_field_coupling_analytic",
    "scalar_potential_coupling",
    "vector_potential_coupling",
]

import time
import numpy as np


from .utils import get_quad_points
from .mesh_calculus import gradient_matrix, mass_matrix
from .integrals import triangle_potential_dipole_linear, triangle_potential_uniform
from .integrals import (
    triangle_potential_approx,
    potential_vertex_dipoles,
)


def magnetic_field_coupling(mesh, r, Nchunks=None, quad_degree=1, analytic=False):
    """
    Computes the coupling matrix which gives the magnetic field at
    target points due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evalution points
    quad_degree: int >= 1
        Quadrature degree (Dunavant scheme) to use.
    analytic: boolean
        compute field using analytic formula (True) or quadrature (False)

    Returns
    -------
    C: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to B-field at the evaluation points

    """

    if analytic:
        return magnetic_field_coupling_analytic(mesh, r, Nchunks)

    mu0 = 4 * np.pi * 1e-7
    coef = mu0 / (4 * np.pi)

    print(
        "Computing magnetic field coupling matrix, %d vertices by %d target points... "
        % (len(mesh.vertices), len(r)),
        end="",
    )
    start = time.time()

    w_quad, r_quad = get_quad_points(
        mesh.vertices, mesh.faces, method="dunavant_0" + str(quad_degree)
    )

    # Rotated gradients (currents)
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)

    # Initialize C-matrix
    n_target_points = len(r)
    n_verts = len(mesh.vertices)
    C = np.zeros((n_target_points, n_verts, 3))

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0] // 100
        else:
            Nchunks = 1

    for n in range(Nchunks):
        # Diffence vectors (Neval, Ntri, Nquad, 3)
        RR = r_quad[None, :, :, :] - r[n::Nchunks, None, None, :]

        # RR/norm(RR)**3 "Gradient of Green's function"
        g = -RR / ((np.linalg.norm(RR, axis=-1) ** 3)[:, :, :, None])

        # Sum over quad points and multiply by triangle area
        g = (g * w_quad[:, None]).sum(axis=-2)
        g *= mesh.area_faces[:, None]

        # Cross product RR/norm(RR)
        C[n::Nchunks, :, 0] = g[:, :, 2] @ Gy - g[:, :, 1] @ Gz
        C[n::Nchunks, :, 1] = g[:, :, 0] @ Gz - g[:, :, 2] @ Gx
        C[n::Nchunks, :, 2] = g[:, :, 1] @ Gx - g[:, :, 0] @ Gy

    duration = time.time() - start
    print("took %.2f seconds." % duration)

    C *= coef
    #    return np.moveaxis(C, 2, 1)
    return np.swapaxes(C, 2, 1)


def magnetic_field_coupling_analytic(mesh, r, Nchunks=None):
    """
    Computes the coupling matrix which gives the magnetic field at
    target points due to currents (stream function) on a surface mesh using analytical formulas.

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points
    Nchunks: int
        number of chunks used in the calculation for saving memory

    Returns
    -------
    C: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to B-field at the evaluation points

    """
    from .integrals import omega, gamma0

    coef = 1e-7

    print(
        "Computing magnetic field coupling matrix analytically, %d vertices by %d target points... "
        % (len(mesh.vertices), len(r)),
        end="",
    )
    start = time.time()

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0] // 100
        else:
            Nchunks = 1

    tn = mesh.face_normals

    # Nfaces, 3, 3
    rfaces = mesh.vertices[mesh.faces]

    # Calculate potentials and related coefficients
    gamma_terms = np.zeros((r.shape[0], mesh.faces.shape[0], 3))
    omega_terms = np.zeros((r.shape[0], mesh.faces.shape[0]))
    # Edges Nfaces, 3, 3
    edges = np.roll(rfaces, 1, -2) - np.roll(rfaces, 2, -2)
    for n in range(Nchunks):
        RRchunk = r[n::Nchunks, None, None, :] - rfaces[None, :, :, :]
        # Neval, Nfaces, xyz
        gamma_terms[n::Nchunks] = -np.einsum("nfe,fei->nfi", gamma0(RRchunk), edges)
        omega_terms[n::Nchunks] = omega(RRchunk)

    # 3 (Nfaces, Nverts) sparse matrices
    G = gradient_matrix(mesh, rotated=False)
    R = gradient_matrix(mesh, rotated=True)
    C = np.zeros((3, r.shape[0], mesh.vertices.shape[0]))

    # Accumulate elements by sparse matrix products for x,y, and z components
    for fcomp in range(3):
        C[fcomp] = omega_terms @ G[fcomp]
        # Accumulate gamma terms for each vertex in the triangle
        for gcomp in range(3):
            # Edge @ Rotated_gradient "==" c_coeff
            # Multiplying with R-matrices takes care of c_coeff calculation
            # and accumulation to right vertex
            C[fcomp] += (tn[:, fcomp] * gamma_terms[:, :, gcomp]) @ R[gcomp]

    duration = time.time() - start
    print("took %.2f seconds." % duration)

    C *= coef
    #    return np.moveaxis(C, 0, 1)
    return np.swapaxes(C, 0, 1)


def scalar_potential_coupling(
    mesh, r, Nchunks=None, multiply_coeff=False, approx_far=False, margin=3
):
    """
    Coupling matrix from a stream function on a mesh
    to scalar potential using analytic integrals.

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points
    Nchunks: int
        number of chunks used in the calculation for saving memory
    multiply_coeff: boolean
        If True, multiply result by mu_0/(4*pi)
    approx_far: boolean,
        approximate the potential using simple quadrature
        (see integrals.potential_dipoles) for points far from the source triangles
    margin: float
        cut-off distance for "far" points measured in mean triangle side length.

    Returns
    -------
    U: ndarray (Np, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to scalar potential at the evaluation points
    """

    print(
        "Computing scalar potential coupling matrix, %d vertices by %d target points... "
        % (len(mesh.vertices), len(r)),
        end="",
    )
    start = time.time()

    # Source and eval locations
    R1 = mesh.vertices[mesh.faces]
    R2 = r

    R2chunks, ichunks = get_chunks(R2, Nchunks, True)
    Uf = np.zeros((R2.shape[0], mesh.faces.shape[0], 3))
    far_chunks = []

    for ichunk, R2chunk in zip(ichunks, R2chunks):
        RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
        RRchunk_verts = R2chunk[:, None, :] - mesh.vertices[None, :, :]
        if approx_far:
            temp = np.zeros(RRchunk.shape[:3])
            # near, far = _split_by_distance(mesh, RRchunk, margin)
            near_v, far_v = _split_by_distance(mesh, RRchunk_verts, margin)
            near = mesh.faces_sparse.T @ near_v
            far_chunks.append(far_v)
            # far = np.invert(near)
            temp[:, near, :] = triangle_potential_dipole_linear(
                RRchunk[:, near], mesh.face_normals[near], mesh.area_faces[near]
            )
            # This far approximation does not speed up the computation much
            # because the quadrature points are so many
            # temp[:,far,:] = potential_dipoles(RRchunk[:, far],
            # mesh.face_normals[far],
            # mesh.area_faces[far])
            Uf[ichunk] = temp
        else:
            Uf[ichunk] = triangle_potential_dipole_linear(
                RRchunk, mesh.face_normals, mesh.area_faces
            )

    # Sparse products are equivalent to this
    # Uv = np.zeros((R2.shape[0], mesh.vertices.shape[0]))
    # for ind_f, f in enumerate(mesh.faces):
    #     Uv[:, f] += Uf[:, ind_f]
    from scipy.sparse import csc_matrix

    Nf = len(mesh.faces)
    Nv = len(mesh.vertices)
    M0 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 0])), (Nf, Nv))
    M1 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 1])), (Nf, Nv))
    M2 = csc_matrix((np.ones(Nf), (np.arange(Nf), mesh.faces[:, 2])), (Nf, Nv))
    Uv = Uf[:, :, 0] @ M0 + Uf[:, :, 1] @ M1 + Uf[:, :, 2] @ M2

    # Calcuate far points by vertex based approximation
    if approx_far:
        areas = mass_matrix(mesh, lumped=True).diagonal()
        for ichunk, R2chunk, far in zip(ichunks, R2chunks, far_chunks):
            RRchunk_verts = R2chunk[:, None, :] - mesh.vertices[None, far, :]
            mask = ichunk[:, None] * far
            Uv[mask] = potential_vertex_dipoles(
                RRchunk_verts, mesh.vertex_normals[far], areas[far]
            ).ravel()

    duration = time.time() - start
    print("took %.2f seconds." % duration)

    if multiply_coeff:
        coeff = 1e-7  # mu_0/(4*pi)
    else:
        coeff = 1 / (4 * np.pi)
    return Uv * coeff


def _triangle_coupling(
    mesh, r, Nchunks=None, approx_far=True, margin=2, chunk_clusters=False, planar=False
):
    """

    Parameters
    ----------
    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points
    Nchunks: int
        number of chunks used in the calculation for saving memory
    approx_far : boolean
        speed up the calculation by approxmating far points. The default is True.
    margin : boolean, optional
        defintion of far points in average triangle sidelength. The default is 2.
    chunk_clusters : boolean, optional
        make chunks clusters, may speed up the calculation. The default is False.
    planar : boolean, optional
        Prapagated to triangle_potential_uniform. For planar meshes
        the calculation can be speeded up. The default is False.

    Returns
    -------
    M : TYPE
        DESCRIPTION.

    """

    # Source and eval locations
    R1 = mesh.vertices[mesh.faces]
    R2 = r

    R2chunks, ichunks = get_chunks(R2, Nchunks, chunk_clusters)

    M = np.zeros((R2.shape[0], mesh.faces.shape[0]))
    print("Computing triangle-coupling matrix")

    if planar:
        print("Assuming the mesh is planar (if not, set planar=False)")

    for ichunk, R2chunk in zip(ichunks, R2chunks):
        RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
        if approx_far:
            RRchunk_centers = R2chunk[:, None, :] - mesh.triangles_center[None, :, :]
            temp = np.zeros(RRchunk.shape[:2])
            near, far = _split_by_distance(mesh, RRchunk_centers, margin)
            temp[:, near] = triangle_potential_uniform(
                RRchunk[:, near], mesh.face_normals[near], planar
            )
            temp[:, far] = triangle_potential_approx(
                RRchunk_centers[:, far], mesh.area_faces[far], reg=0
            )
            M[ichunk] = temp
        else:
            M[ichunk] = triangle_potential_uniform(RRchunk, mesh.face_normals, planar)

    return M


def vector_potential_coupling(
    mesh, r, Nchunks=None, approx_far=True, margin=2, chunk_clusters=False
):
    """
    Compute vector potential coupling matrices
    from a linear stream function on a mesh using analytic integrals.

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points
    approx_far: Boolean (True)
        If True, use approximate calculation for triangles that
        far from the source triangles using a simple quadrature
        (see integrals.triangle_potential_approx)
    margin: float
        Cut-off distance for "far" points measured in mean triangle side length.

    Returns
    -------
    A: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to vector potential at the evaluation points
    """

    coeff = 1e-7  # mu_0/(4*pi)

    Af = _triangle_coupling(mesh, r, Nchunks, approx_far, margin, chunk_clusters)

    # Rotated gradients (currents)
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
    # Accumulate the elements
    Av = np.array([Af @ Gx, Af @ Gy, Af @ Gz])

    return Av * coeff


def get_chunks(r, Nchunks, clusters=True):
    """ Chunk points in 'r' to Nchunks

        r : ndarray (Npoints, 3)
    """
    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0] // 100
        else:
            Nchunks = 1
    if clusters:
        # Voronoi cells of random vertices
        i_samples = np.random.randint(0, r.shape[0], Nchunks)
        dists = np.linalg.norm(r[:, None, :] - r[None, i_samples, :], axis=-1)
        labels = np.argmin(dists, axis=1)
        # indices as boolean arrays
        # Number of unique labels can be smaller than Nchunks if
        # there are vertices without any points in their cells
        ichunks = [labels == label for label in np.unique(labels)]
        rchunks = [r[mask] for mask in ichunks]
    else:
        # Chunk r by array split and get the corresponding indices as slices
        rchunks = np.array_split(r, Nchunks, axis=0)
        lengths = [len(ri) for ri in rchunks]
        inds = np.cumsum([0] + lengths)
        ichunks = [slice(inds[i], inds[i + 1]) for i in range(len(lengths))]

    return rchunks, ichunks


def _split_by_distance(mesh, RR, margin=3):
    avg_sidelength = np.sqrt(4 / np.sqrt(3) * np.mean(mesh.area_faces[::100]))

    RRnorm = np.linalg.norm(RR, axis=-1)
    # near = np.nonzero(np.min(RRnorm, axis=0) < avg_sidelength * margin)[0]
    # far = np.setdiff1d(np.arange(0, len(mesh.faces)), near, assume_unique=True)
    near = np.min(RRnorm, axis=0) < avg_sidelength * margin
    far = np.invert(near)

    return near, far


# def _split_by_distance(mesh, RR, margin=3):
#     avg_sidelength = np.sqrt(4/np.sqrt(3)*np.mean(mesh.area_faces[::100]))
# #    np.mean(np.linalg.norm(np.diff(mesh.vertices[mesh.edges[::1000]], axis=1), axis=-1))

#     RRnorm = np.linalg.norm(RR, axis=-1)
#     near = np.nonzero(np.min(RRnorm, axis=(0, 2)) < avg_sidelength * margin)[0]

#     far = np.setdiff1d(np.arange(0, len(mesh.faces)), near, assume_unique=True)

# #    print('near: %d, far: %d'%(len(near), len(far)))
#     return near, far
