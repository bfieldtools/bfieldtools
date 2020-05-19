"""
Contains functions for computing thermal noise in conductive thin objects.

"""

__all__ = [
    "compute_current_modes",
    "noise_covar",
    "noise_covar_dir",
    "noise_var",
    "visualize_current_modes",
]

import numpy as np
import trimesh

from .suhtools import SuhBasis
from .mesh_conductor import MeshConductor


def compute_current_modes(
    obj,
    T,
    resistivity,
    thickness,
    mode="AC",
    freqs=np.array((0,)),
    Nmodes=None,
    return_eigenvals=False,
    **kwargs
):
    """
    Calculates the (AC or DC) Johnson noise current modes on the conducting surface.
    
    Parameters
    ----------
    obj: Trimesh-object or MeshConductor-object
        Represents the boundary on which current density is specified
        or MeshConductor object that wraps the mesh
    T: float
        Temperature in Kelvins
    resistivity: float or array (Nfaces)
        Resistivity value in Ohm/meter
    thickness: float or array (Nfaces)
        Thickness of the surface. NB! Must be small in comparison to observation distance
    mode: 'AC' or 'DC'
        Calculate modes as a function of frequency or just at DC?
    freqs: Nfreqs array
        The frequencies at which the eddy-current modes are computed. Obsolete for 
        'DC' mode. In 'AC', calculate at DC in default.
    Nmodes: int
        How many modes are computed? If None, all Nvertices modes are computed
    return_eigenvals: boolean
        Return also the eigenvalues (the inverse circuit time constants)?
    kwargs: dict
        Passed to Conductor creation if a Trimesh object is passed as 'obj'
        
    Returns
    -------
    vl: (Nvertices x Nmodes x Nfreqs) (AC) or (Nvertices x Nmodes) (DC) array
        The spectral eddy-current modes
    u: Nmodes array
        The eigenvalues
    """

    kB = 1.38064852e-23

    if isinstance(obj, MeshConductor):
        obj.resistivity = resistivity
        obj.thickness = thickness
        print(
            "Updated MeshConductor object thickness and resistivity attributes according to function parameters"
        )
    elif isinstance(obj, trimesh.Trimesh):
        obj = MeshConductor(
            mesh_obj=obj,
            resistivity=resistivity,
            thickness=thickness,
            resistance_full_rank=False,
            **kwargs
        )
    else:
        raise TypeError("obj type should be either Trimesh or Conductor")

    if mode == "AC":
        suh = SuhBasis(obj, Nc=Nmodes, magnetic="AC")
    elif mode == "DC":
        suh = SuhBasis(obj, Nc=Nmodes, magnetic="DC")
    else:
        raise ValueError("Mode should be either 'AC' or 'DC'")

    v = suh.basis
    u = suh.eigenvals

    Nfreqs = len(freqs)

    if mode == "AC":
        Nfreqs = len(freqs)
        vl = np.zeros((len(suh.mesh_conductor.mesh.vertices), v.shape[1], Nfreqs))

        # Scale the eigenmodes with the spectral density of the thermal noise current
        for i in range(v.shape[1]):
            amp = (
                2
                * np.sqrt(kB * T / u[i])
                * np.sqrt(1 / (1 + (2 * np.pi * freqs / u[i]) ** 2))
            )
            vl[:, i, :] = (
                (
                    np.zeros((Nfreqs, vl.shape[0]))
                    + suh.mesh_conductor.inner2vert @ v[:, i]
                ).T
            ) * amp
    elif mode == "DC":
        vl = np.zeros((len(suh.mesh_conductor.mesh.vertices), v.shape[1]))

        # Scale the eigenmodes with the spectral density of the thermal noise current
        for i in range(v.shape[1]):
            amp = 2 * np.sqrt(kB * T / u[i])
            vl[:, i] = suh.mesh_conductor.inner2vert @ v[:, i] * amp

    if return_eigenvals:
        return vl, u
    else:
        return vl


def noise_covar(B_coupling, vl, Nmodes=None):
    """
    Calculates (AC or DC) magnetic noise covariance along x, y and z from the 
    modes vl.
    
    Parameters
    ----------
    B_coupling: ndarray (Np, 3, Nvertices)
        Magnetic field coupling matrix from the mesh
    vl: ndarray (Nvertices, Nmodes, x Nfreqs) or (Nvertices, Nmodes)
        The Johnson noise current modes on the mesh
    Nmodes: int
        How many modes are included? If None, all modes in vl are included

    Returns
    -------
    Bcov: ndarray (Np, Np, 3, Nfreqs) or (Np, Np, 3)
        Magnetic noise covariance
    """
    if Nmodes is None:
        Nmodes = vl.shape[1]

    if vl.ndim == 2:
        b = np.einsum("ihj,jl->ilh", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("jih,lih->jlh", b, b)
    else:
        b = np.einsum("ihj,jlk->ilhk", B_coupling, vl[:, 0:Nmodes, :])
        Bcov = np.einsum("jihk,lihk->jlhk", b, b)

    return Bcov


def noise_var(B_coupling, vl, Nmodes=None):
    """
    Calculates (AC or DC) magnetic noise variance along x, y and z from the 
    modes vl.
    
    Parameters
    ----------
    B_coupling: ndarray (Np, 3, Nvertices)
        Magnetic field coupling matrix from the mesh
    vl: ndarray (Nvertices, Nmodes, x Nfreqs) or (Nvertices, Nmodes)
        The Johnson noise current modes on the mesh
    Nmodes: int
        How many modes are included? If None, all modes in vl are included

    Returns
    -------
    Bcov: ndarray (Np, 3, Nfreqs) or (Np, 3)
        Magnetic noise variance
    """
    if Nmodes is None:
        Nmodes = vl.shape[1]

    if vl.ndim == 2:
        b = np.einsum("ihj,jl->ilh", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("ijh,ijh->ih", b, b)
    else:
        b = np.einsum("ihj,jlk->ilhk", B_coupling, vl[:, 0:Nmodes, :])
        Bcov = np.einsum("ijhk,ijhk->ihk", b, b)

    return Bcov


def noise_covar_dir(B_coupling, vl, Nmodes=None):
    """
    Calculates (AC or DC) magnetic noise covariance between x, y and z directions
    from the modes vl.
    
    Parameters
    ----------
    B_coupling: ndarray (Np, 3, Nvertices)
        Magnetic field coupling matrix from the mesh
    vl: ndarray (Nvertices, Nmodes, x Nfreqs) or (Nvertices, Nmodes)
        The Johnson noise current modes on the mesh
    Nmodes: int
        How many modes are included? If None, all modes in vl are included

    Returns
    -------
    Bcov: ndarray (Np, 3, 3, Nfreqs) or (Np, 3, 3)
        Magnetic noise covariance x, y and z field components
    """
    if Nmodes is None:
        Nmodes = vl.shape[1]

    if vl.ndim == 2:
        b = np.einsum("ihj,jl->ilh", B_coupling, vl[:, 0:Nmodes])
        Bcov = np.einsum("ihj,ihl-> ijl", b, b)
    else:
        b = np.einsum("ihj,jlk->ilhk", B_coupling, vl[:, 0:Nmodes, :])
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
    from mayavi import mlab

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
