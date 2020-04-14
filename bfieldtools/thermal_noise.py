"""
Contains functions for computing thermal noise in conductive thin objects.

"""

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh  # , eigsh
from mayavi import mlab

from .mesh_magnetics import magnetic_field_coupling
from .mesh_calculus import laplacian_matrix, mass_matrix
from .mesh_properties import self_inductance_matrix, resistance_matrix
from . import utils


def compute_AC_current_modes(
    mesh, M, R, freqs, T, closed=True, Nmodes=None, return_eigenvals=False
):
    """
    Parameters
    ----------
    mesh: Trimesh mesh object
        The surface mesh
    M: (Nvertices x Nvertices) array
        The self-inductance matrix of `mesh`
    R: (Nvertices x Nvertices) array
        The resistance matrix of `mesh`
    freqs: Nfreqs array
        The frequencies at which the eddy-current modes are computed
    T: float
        Temperature in Kelvins
    closed: boolean
        Is the mesh closed (True) or not (False)
    Nmodes: int
        How many modes are computed? If None, all Nvertices modes are computed
    return_eigenvals: boolean
        Return also the eigenvalues (the inverse circuit time constants)?

    Returns
    -------
    vl: (Nvertices x Nmodes x Nfreqs) array
        The spectral eddy-current modes
    u: Nmodes array
        The eigenvalues (the inverse circuit time constants)

    """

    kB = 1.38064852e-23

    boundary_verts, inner_verts = utils.find_mesh_boundaries(mesh)

    # Ensure that M is symmetric
    M = 0.5 * (M + M.T)

    R = R.toarray()  # convert R to array

    # If mesh is closed, add 'deflation' to the matrices
    if closed:
        R += np.ones_like(R) * np.mean(np.diag(R))
        M += np.ones_like(M) * np.mean(np.diag(M))

    # Compute the eigenmodes
    if Nmodes == None:
        u, v = eigh(R[inner_verts][:, inner_verts], M[inner_verts][:, inner_verts])
    else:
        u, v = eigh(
            R[inner_verts][:, inner_verts],
            M[inner_verts][:, inner_verts],
            eigvals=(0, Nmodes),
        )

    Nfreqs = freqs.shape[0]

    # Scale the eigenmodes with the spectral density of the thermal noise current
    vl = np.zeros((M.shape[0], v.shape[1], Nfreqs))

    for i in range(v.shape[1]):
        amp = (
            2
            * np.sqrt(kB * T / u[i])
            * np.sqrt(1 / (1 + (2 * np.pi * freqs / u[i]) ** 2))
        )
        vl[inner_verts, i, :] = ((np.zeros((Nfreqs, v.shape[0])) + v[:, i]).T) * amp

    if return_eigenvals:
        return vl, u
    else:
        return vl


def compute_DC_current_modes(
    mesh, R, T, closed=True, Nmodes=None, return_eigenvals=False
):
    """
    Parameters
    ----------
    mesh: Trimesh mesh object
        The surface mesh
    R: (Nvertices x Nvertices) array
        The resistance matrix of `mesh`
    T: float
        Temperature in Kelvins
    closed: boolean
        Is the mesh closed (True) or not (False)
    Nmodes: int
        How many modes are computed? If None, all Nvertices modes are computed
    return_eigenvals: boolean
        Return also the eigenvalues (the inverse circuit time constants)?

    Returns
    -------
    vl: (Nvertices x Nmodes x Nfreqs) array
        The spectral eddy-current modes
    u: Nmodes array
        The eigenvalues (the inverse circuit time constants)

    """

    kB = 1.38064852e-23

    boundary_verts, inner_verts = utils.find_mesh_boundaries(mesh)

    R = R.toarray()  # convert R to array

    # If mesh is closed, add 'deflation' to the matrices
    if closed:
        R += np.ones_like(R) * np.mean(np.diag(R))

    M = mass_matrix(mesh)
    # Compute the eigenmodes
    if Nmodes == None:
        u, v = eigh(
            R[inner_verts][:, inner_verts], M.todense()[inner_verts][:, inner_verts]
        )
    else:
        u, v = eigh(
            R[inner_verts][:, inner_verts],
            M.todense()[inner_verts][:, inner_verts],
            eigvals=(0, Nmodes),
        )

    # Scale the eigenmodes with the spectral density of the thermal noise current
    vl = np.zeros((M.shape[0], v.shape[1]))

    for i in range(v.shape[1]):
        amp = 2 * np.sqrt(kB * T / u[i])
        vl[inner_verts, i] = v[:, i] * amp

    if return_eigenvals:
        return vl, u
    else:
        return vl


def compute_current_modes(mesh, boundaries=None, return_eigenvals=False):
    """
    Computes eddy-current modes for a mesh using surface laplacian.
    Uses Dirichlet boundary condition, i.e., stream function is zero at boundary:
    no current flow outside the surface.
    The modes are normalized so that the squared norm of the stream function gradient
    integrates to 1 over the surface. With this normalization, the resistances
    of the current modes are R_k = 1/(sigma*d), sigma = conductivity, d = thickness.
    See Zevenhoven et al. (2014).

    Parameters
    ----------
    mesh: Trimesh mesh object
        The surface mesh
    boundaries: list of N_holes

    Returns
    -------
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]

    """
    boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(
        mesh.vertices, mesh.faces, mesh.edges
    )

    if boundaries:
        L_holes = laplacian_matrix(mesh, None, inner_verts, boundaries)
        M_holes = mass_matrix(mesh, inner_verts, boundaries)

        u, v = eigh(-L_holes.todense(), M_holes.todense())

        # Normalize the laplacien eigenvectors

        for i in range(v.shape[1]):
            v[:, i] = v[:, i] / np.sqrt(u[i])

        # Assign values per vertex
        vl = np.zeros((mesh.vertices.shape[0], v.shape[1]))

        vl[inner_verts] = v[: -len(boundaries)]

        for b_idx, b in enumerate(boundaries):
            vl[b] = v[len(inner_verts) + b_idx]

    else:
        L = laplacian_matrix(mesh)
        M = mass_matrix(mesh)

        u, v = eigh(
            -L.todense()[inner_verts][:, inner_verts],
            M.todense()[inner_verts][:, inner_verts],
        )

        # Normalize the laplacien eigenvectors
        vl = np.zeros(M.shape)
        for i in range(v.shape[1]):
            vl[inner_verts, i] = v[:, i] / np.sqrt(u[i])

    if return_eigenvals:
        return vl, u
    else:
        return vl


def noise_covar(mesh, B_coupling, vl, Nmodes=None):
    if Nmodes == None:
        Nmodes = vl.shape[1]

    if vl.ndim == 2:
        b = np.einsum("ihj,jl->ilh", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("jih,lih->jlh", b, b)
    else:
        b = np.einsum("ihj,jlk->ilhk", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("jihk,lihk->jlhk", b, b)

    return Bcov


def noise_var(mesh, B_coupling, vl, Nmodes=None):
    if Nmodes == None:
        Nmodes = vl.shape[1]

    if vl.ndim == 2:
        b = np.einsum("ihj,jl->ilh", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("ijh,ijh->ih", b, b)
    else:
        b = np.einsum("ihj,jlk->ilhk", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("ijhk,ijhk->ihk", b, b)

    return Bcov


def noise_covar_dir(mesh, B_coupling, vl, Nmodes=None):
    if Nmodes == None:
        Nmodes = vl.shape[1]

    if vl.ndim == 2:
        b = np.einsum("ihj,jl->ilh", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("ihj,ihl-> ijl", b, b)
    else:
        b = np.einsum("ihj,jlk->ilhk", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("ihjk,ihlk-> ijlk", b, b)

    return Bcov


def visualize_current_modes(
    mesh, vl, Nmodes, scale, contours=True, colormap="bwr", dist=0.5
):
    """
    Visualizes current modes up to Nmodes.

    Parameters
    ----------
    mesh: Trimesh mesh object
        The surface mesh
    vl: Nvertices x Nvertices array
        The normalized eddy-current modes vl[:,i]
    Nmodes: int
        Number of modes to be plotted
    scale: float
        Scaling factor
    contours: boolean
        If True, show contours
    colormap: string
        Which (matplotlib) colormap to use

    """

    N1 = np.floor(np.sqrt(Nmodes))
    dx = (mesh.vertices[:, 0].max() - mesh.vertices[:, 0].min()) * (1 + dist)
    dy = (mesh.vertices[:, 1].max() - mesh.vertices[:, 1].min()) * (1 + dist)

    i = 0
    j = 0
    for n in range(Nmodes):
        print(i, j)
        points = mesh.vertices.copy()
        points[:, 0] += i * dx
        points[:, 1] += j * dy
        s = mlab.triangular_mesh(
            *points.T, mesh.faces, scalars=vl[:, n], colormap=colormap
        )

        limit = np.max(np.abs(vl[:, n]))

        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit, limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.enable_contours = contours

        if i < N1:
            i += 1
        else:
            j += 1
            i = 0

    return s
