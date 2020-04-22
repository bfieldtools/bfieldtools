import numpy as np
from mayavi import mlab

from bfieldtools import suhtools
from bfieldtools.utils import load_example_mesh

from .test_conductor import _fake_conductor

import pytest


def test_suhbasis():

    mesh = load_example_mesh("unit_disc")

    for obj in [mesh, _fake_conductor(mesh_name="unit_disc")]:
        for mag in [False]:
            for bc in ["neumann", "dirichlet"]:
                suh = suhtools.SuhBasis(obj, Nc=10, boundary_condition=bc, magnetic=mag)

                suh.field(suh.basis[1], np.array([[0, 0, 1]]))

                suh.plot(Nfuncs=3)

    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="dirichlet", magnetic="DC")
    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="dirichlet", magnetic="AC")

    try:
        suh = suhtools.SuhBasis(10, Nc=10)
    except:
        print("Caught test exception")

    f = mlab.figure()

    suh.plot([0, 2], figure=f)
    suh.plot(Nfuncs=4, Ncols=2, figure=f)


@pytest.mark.xfail
def test_suhbasis_fail():
    mesh = load_example_mesh("unit_disc")
    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="foo")
