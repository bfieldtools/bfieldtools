import numpy as np
import os
from bfieldtools import kicad_io
from bfieldtools.utils import load_example_mesh

from .test_mesh_conductor import _fake_mesh_conductor, _fake_streamfunction

import pytest


def test_kicad_io():

    s = _fake_streamfunction()
    c = s.mesh_conductor

    lp = s.discretize(6)

    kicad_io.python_to_kicad(
        [loop.discrete(lp.vertices) for loop in lp.entities],
        "test_file",
        plane_axes=(0, 1),
        origin=(-1, -1),
        layer="F.Cu",
        net=1,
        scaling=1,
    )

    os.remove("test_file")
