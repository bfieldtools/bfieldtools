import numpy as np

from bfieldtools import thermal_noise
from bfieldtools.utils import load_example_mesh

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
