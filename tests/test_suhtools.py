import numpy as np

from bfieldtools import suhtools
from bfieldtools.utils import load_example_mesh

from .test_conductor import _fake_conductor

import pytest


def test_suhbasis():

    mesh = load_example_mesh("unit_disc")

    for obj in [mesh, _fake_conductor(mesh_name="unit_disc")]:
        for mag in [False]:
            for bc in ["neumann", "dirichlet"]:
                suh = suhtools.SuhBasis(
                    mesh, Nc=10, boundary_condition=bc, magnetic=mag
                )

                suh.field(suh.basis[1], np.array([[0, 0, 1]]))

    suh = suhtools.SuhBasis(mesh, Nc=10, boundary_condition="dirichlet", magnetic=True)

    try:
        suh = suhtools.SuhBasis(10, Nc=10)
    except:
        print("Caught test exception")

    suh.plot(Nfuncs=3)
