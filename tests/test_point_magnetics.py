# -*- coding: utf-8 -*-
"""
Test point magnetics

"""

from bfieldtools.point_magnetics import *


def test_magnetic_field():
    p = np.array(
        [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]]
    )

    p0 = points / 2
    m = np.array([0, 0, 1]) * np.ones((5, 1))

    field = magnetic_field_dipoles(p, m, p0)

    # TODO analytic solutions


def test_vector_potential():
    p = np.array(
        [[-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [-1.0, 1.0, 0.0]]
    )

    p0 = points / 2
    m = np.array([0, 0, 1]) * np.ones((5, 1))
    field = vector_potential_dipoles(p, m, p0)
