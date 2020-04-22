"""
Contains functions for computing thermal noise in conductive thin objects.

"""

import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh  # , eigsh
from mayavi import mlab

from .suhtools import SuhBasis
from .mesh_magnetics import magnetic_field_coupling
from .mesh_calculus import laplacian_matrix, mass_matrix
from .mesh_impedance import self_inductance_matrix, resistance_matrix
from . import utils


def compute_AC_current_modes(
    obj, freqs, T, closed=True, Nmodes=None, return_eigenvals=False
):
    """
    Parameters
    ----------
    obj : Trimesh-object or Conductor-object
        Represents the boundary on which current density is specified
        or Conductor object that wraps the mesh
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

    suh = SuhBasis(obj, Nc=Nmodes, magnetic="AC")

    v = suh.basis
    u = suh.eigenvals

    Nfreqs = freqs.shape[0]

    # Scale the eigenmodes with the spectral density of the thermal noise current
    vl = np.zeros((len(suh.conductor.mesh.vertices), v.shape[1], Nfreqs))

    for i in range(v.shape[1]):
        amp = (
            2
            * np.sqrt(kB * T / u[i])
            * np.sqrt(1 / (1 + (2 * np.pi * freqs / u[i]) ** 2))
        )
        vl[suh.inner_vertices, i, :] = (
            (np.zeros((Nfreqs, v.shape[0])) + v[:, i]).T
        ) * amp

    if return_eigenvals:
        return vl, u
    else:
        return vl


def compute_DC_current_modes(obj, T, closed=True, Nmodes=None, return_eigenvals=False):
    """
    Parameters
    ----------
    obj : Trimesh-object or Conductor-object
        Represents the boundary on which current density is specified
        or Conductor object that wraps the mesh
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

    suh = SuhBasis(obj, Nc=Nmodes, magnetic="DC")

    v = suh.basis
    u = suh.eigenvals

    # Scale the eigenmodes with the spectral density of the thermal noise current
    vl = np.zeros((len(suh.conductor.mesh.vertices), v.shape[1]))

    for i in range(v.shape[1]):
        amp = 2 * np.sqrt(kB * T / u[i])
        vl[suh.inner_vertices, i] = v[:, i] * amp

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
