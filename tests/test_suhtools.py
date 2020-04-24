import numpy as np
from mayavi import mlab

from bfieldtools import suhtools
from bfieldtools.utils import load_example_mesh

from .test_mesh_conductor import _fake_mesh_conductor

import pytest

from numpy.testing import assert_allclose


def test_suhbasis():

    mesh = load_example_mesh("unit_disc")

    for obj in [mesh, _fake_mesh_conductor(mesh_name="unit_disc")]:
        for mag in [False]:
            for bc in ["neumann", "dirichlet"]:
                suh = suhtools.SuhBasis(obj, Nc=10, boundary_condition=bc, magnetic=mag)
                assert suh.basis.shape[1] == 10

                suh.field(suh.basis[1], np.array([[0, 0, 1]]))

                suh.plot(Nfuncs=3)

    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="dirichlet", magnetic="DC")
    assert suh.basis.shape[1] == 10

    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="dirichlet", magnetic="AC")
    assert suh.basis.shape[1] == 10

    suh = suhtools.SuhBasis(mesh, boundary_condition="dirichlet", magnetic="DC")
    suh = suhtools.SuhBasis(mesh, boundary_condition="dirichlet", magnetic="AC")
    suh = suhtools.SuhBasis(mesh, boundary_condition="dirichlet", magnetic=False)

    suh.calculate_basis(shiftinvert=False, v0=np.ones(57))

    try:
        suh = suhtools.SuhBasis(10, Nc=10)
        assert suh.basis.shape[1] == 10
    except:
        print("Caught test exception")

    f = mlab.figure()

    suh.plot([0, 2], figure=f)
    suh.plot(Nfuncs=4, Ncols=2, figure=f)


def test_shiftinvert():

    mesh = load_example_mesh("unit_disc")

    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="dirichlet", magnetic=False)

    suh.calculate_basis(shiftinvert=True)
    shiftinvert_eigenvals = suh.eigenvals
    shiftinvert_basis = suh.basis

    suh.calculate_basis(shiftinvert=False)
    nonshiftinvert_eigenvals = suh.eigenvals
    nonshiftinvert_basis = suh.basis

    assert_allclose(shiftinvert_eigenvals, nonshiftinvert_eigenvals)


def test_suhbasis_closed():
    import trimesh

    # Low-res ico-sphere, but need to copy since Trimesh primitives are protected
    mesh = trimesh.primitives.Sphere(subdivisions=2)
    mesh = trimesh.Trimesh(mesh.vertices, mesh.faces)

    for mag in [False, "DC", "AC"]:
        suh = suhtools.SuhBasis(mesh, Nc=10, magnetic=mag)
        assert suh.basis.shape[1] == 10

        suh = suhtools.SuhBasis(mesh, magnetic=mag)
        assert suh.basis.shape[1] == suh.Nc

    suh = suhtools.SuhBasis(mesh, Nc=41, magnetic=mag)
    assert suh.basis.shape[1] == suh.Nc


@pytest.mark.xfail
def test_suhbasis_fail():
    mesh = load_example_mesh("unit_disc")
    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="foo")


@pytest.mark.xfail
def test_suhbasis_fail2():
    mesh = load_example_mesh("unit_disc")
    suh = suhtools.SuhBasis(mesh, Nc=10, magnetic="foo")


@pytest.mark.xfail
def test_suhbasis_fail3():
    mesh = load_example_mesh("unit_disc")
    suh = suhtools.SuhBasis(mesh, Nc=1000)
