from bfieldtools import conductor
from bfieldtools.utils import load_example_mesh

import pytest

import numpy as np
from numpy.testing import (
    assert_array_almost_equal,
    assert_array_equal,
    assert_allclose,
    assert_equal,
)


def _fake_conductor(mesh_name="10x10_plane", **kwargs):
    """
    Creates an example 'fake' Conductor object
    """
    return conductor.Conductor(mesh_obj=load_example_mesh(mesh_name), **kwargs)


def _fake_streamfunction(vals=None, mesh_name="unit_disc", **kwargs):
    """
    Creates an example 'fake' StreamFunction object
    """
    if vals:
        return conductor.StreamFunction(vals, _fake_conductor(mesh_name, **kwargs))
    else:
        mesh = load_example_mesh(mesh_name)

        vals = 1 - np.linalg.norm(mesh.vertices, axis=1)

        return conductor.StreamFunction(
            vals, _fake_conductor(mesh_name, basis_name="vertex")
        )


def test_conductor_creation():
    """
    Tests different ways to create Conductor objects, check that basis operators work
    """
    for test_mesh in ["unit_sphere", "unit_disc", "plane_w_holes"]:
        for basis_name in ["suh", "inner", "vertex"]:

            c = _fake_conductor(mesh_name=test_mesh, basis_name=basis_name, N_suh=10)

            assert c.inner2vert.shape == (
                len(c.mesh.vertices),
                len(c.inner_vertices) + len(c.holes),
            )
            assert c.vert2inner.shape == (
                len(c.inner_vertices) + len(c.holes),
                len(c.mesh.vertices),
            )

            assert_array_almost_equal(
                (c.vert2inner @ c.inner2vert).toarray(),
                np.identity(len(c.inner_vertices) + len(c.holes)),
            )

            inner_diag = np.zeros((len(c.mesh.vertices), len(c.mesh.vertices)))
            inner_diag[c.inner_vertices, c.inner_vertices] = 1

            for hole in c.holes:
                inner_diag[np.asarray(hole)[:, None], np.asarray(hole)] += 1 / len(hole)

            assert_array_equal((c.inner2vert @ c.vert2inner).toarray(), inner_diag)

            if basis_name == "suh":
                assert c.basis.shape == (len(c.mesh.vertices), c.opts["N_suh"])


def test_streamfunction_creation():
    """
    Test creating StreamFunctions with different bases
    """

    conductor.StreamFunction(
        np.zeros((10,)), _fake_conductor(basis_name="suh", N_suh=10)
    )
    conductor.StreamFunction(np.zeros((584,)), _fake_conductor(basis_name="inner"))
    conductor.StreamFunction(np.zeros((676,)), _fake_conductor(basis_name="vertex"))


def test_conductor_attributes():
    """
    tests Conductor attributes
    """
    c = _fake_conductor("unit_disc")

    c2 = _fake_conductor("unit_disc")
    c2.mesh.vertices += np.array([0, 0, 1])

    ind = c.inductance
    m_ind = c.mutual_inductance(c2)

    c2.sph_couplings

    # This should cause an error
    try:
        c.set_basis("foo")
    except:
        print("Failed as expected")

    c.plot_mesh()

    c.U_coupling(np.array([[0, 0, 2]]))
    c.U_coupling.reset()


def test_streamfunction_attributes():
    """
    tests StreamFunction functionality
    """
    s = _fake_streamfunction(mesh_name="unit_disc")
    s.__repr__()

    assert len(s.inner) == len(s.conductor.inner_vertices)
    assert len(s.vert) == len(s.conductor.mesh.vertices)

    s.power
    s.magnetic_energy

    s.plot()
    s.plot(contours=6)

    lp = s.discretize(N_contours=6)
