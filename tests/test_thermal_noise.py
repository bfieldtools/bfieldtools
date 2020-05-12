import numpy as np

from bfieldtools import thermal_noise
from bfieldtools.utils import load_example_mesh
from bfieldtools.mesh_magnetics import magnetic_field_coupling

import trimesh
import pytest

from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
    assert_allclose,
    assert_equal,
)

from .test_mesh_conductor import _fake_mesh_conductor


def test_compute_modes():

    obj = _fake_mesh_conductor("unit_disc")

    DC_n = thermal_noise.compute_current_modes(
        obj, T=300, resistivity=1e-8, thickness=1e-3, mode="DC"
    )

    AC_n = thermal_noise.compute_current_modes(
        obj, freqs=np.array((0,)), T=300, resistivity=1e-8, thickness=1e-3, mode="AC"
    )

    points = np.array([[0.1, 0, 0.2], [-0.1, 0.1, -0.2]])

    obj.set_basis("vertex")

    DC_Bn_covar = thermal_noise.noise_covar(obj.B_coupling(points), DC_n)
    AC_Bn_covar = thermal_noise.noise_covar(obj.B_coupling(points), AC_n)

    assert_allclose(DC_Bn_covar, AC_Bn_covar[:, :, :, 0])

    DC_Bn_covar_dir = thermal_noise.noise_covar_dir(obj.B_coupling(points), DC_n)
    AC_Bn_covar_dir = thermal_noise.noise_covar_dir(obj.B_coupling(points), AC_n)

    assert_allclose(DC_Bn_covar_dir, AC_Bn_covar_dir[:, :, :, 0])

    DC_Bn_var = thermal_noise.noise_var(obj.B_coupling(points), DC_n)
    AC_Bn_var = thermal_noise.noise_var(obj.B_coupling(points), AC_n)

    assert_allclose(DC_Bn_var, AC_Bn_var[:, :, 0])

    thermal_noise.visualize_current_modes(obj.mesh, DC_n, scale=0.5, Nmodes=3)


def test_vs_analytic():

    obj = _fake_mesh_conductor("unit_disc")
    mesh = obj.mesh

    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    # Fix the simulation parameters
    d = 100e-6  # thickness
    sigma = 3.7e7  # conductivity
    res = 1 / sigma  # resistivity
    T = 300  # temperature
    kB = 1.38064852e-23  # Boltz
    mu0 = 4 * np.pi * 1e-7  # permeability of freespace

    # Compute the AC-current modes and visualize them
    vl, u = thermal_noise.compute_current_modes(
        obj=mesh, T=T, resistivity=res, thickness=d, mode="AC", return_eigenvals=True
    )

    # Define field points on z axis
    Np = 10
    z = np.linspace(0.2, 1, Np)
    fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T

    B_coupling = magnetic_field_coupling(
        mesh, fp, analytic=True
    )  # field coupling matrix

    # Compute noise variance
    B = np.sqrt(thermal_noise.noise_var(B_coupling, vl))

    # Next, we compute the DC noise without reference to the inductance
    vl_dc, u_dc = thermal_noise.compute_current_modes(
        obj=mesh, T=T, resistivity=res, thickness=d, mode="DC", return_eigenvals=True
    )

    # Compute noise variance
    B_dc = np.sqrt(thermal_noise.noise_var(B_coupling, vl_dc))

    # Calculate Bz noise using analytical formula and plot the results
    r = 1
    Ban = (
        mu0
        * np.sqrt(sigma * d * kB * T / (8 * np.pi * z ** 2))
        * (1 / (1 + z ** 2 / r ** 2))
    )

    assert_allclose(B_dc[:, 2], B[:, 2, 0])

    assert_allclose(Ban, B[:, 2, 0], rtol=5e-2)
