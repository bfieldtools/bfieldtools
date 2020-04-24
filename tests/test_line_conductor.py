# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 11:32:38 2020

@author: Rasmus Zetter
"""


from bfieldtools import line_conductor
from bfieldtools.utils import load_example_mesh


import pytest

import numpy as np


def _fake_line_conductor():
    """
    return synthetic LineConductor
    """
    loops = [
        np.array([[0, 0, 1], [0, 1, 0], [0, 2, 0], [1, 2, 1], [1, 1, 0]]),
        2 * np.array([[0, 0, 1], [2, 1, 0], [1, 2, 0], [1, 2, 1], [1, 1, 0]]),
    ]
    return line_conductor.LineConductor(loops)


def _fake_line_conductor_closed():
    """
    return synthetic LineConductor
    """
    loops = [
        np.array([[0, 0, 1], [0, 1, 0], [0, 2, 0], [1, 2, 1], [1, 1, 0], [0, 0, 1]]),
        2
        * np.array([[0, 0, 1], [2, 1, 0], [1, 2, 0], [1, 2, 1], [1, 1, 0], [0, 0, 1]]),
    ]
    return line_conductor.LineConductor(loops)


def test_line_conductor_functionality():
    """
    Test LineConductor creation and functionality
    """

    lp = _fake_line_conductor()
    lp2 = _fake_line_conductor_closed()

    lp.plot_loops()

    mesh = load_example_mesh("unit_disc")
    scalars = np.ones((len(mesh.vertices),)) - np.linalg.norm(mesh.vertices, axis=1)

    lp = line_conductor.LineConductor(mesh=mesh, scalars=scalars, N_contours=6)

    slp = lp.simplify()

    B = slp.magnetic_field(mesh.vertices + np.array([0, 0, 1]))
    A = slp.vector_potential(mesh.vertices + np.array([0, 0, 1]))
    U = slp.scalar_potential(mesh.vertices + np.array([0, 0, 1]))

    mesh.vertices += np.array([[0, 0, 1]])

    M = slp.mesh_mutual_inductance(mesh)
