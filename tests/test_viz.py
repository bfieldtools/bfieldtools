# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:09:56 2020

@author: Rasmus Zetter
"""

import pytest
from bfieldtools import viz
from bfieldtools.utils import load_example_mesh

from .test_conductor import _fake_conductor, _fake_streamfunction


import numpy as np


def test_mesh_plotting():

    mesh = load_example_mesh("unit_disc")

    v_scalars = np.random.rand(len(mesh.vertices))
    f_scalars = np.random.rand(len(mesh.faces))

    viz.plot_mesh(mesh)

    viz.plot_data_on_faces(mesh, f_scalars)
    viz.plot_data_on_vertices(mesh, v_scalars)


def test_plot_cross_section():

    X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))

    Z = 3 - (X ** 2 + Y ** 2)

    viz.plot_cross_section(X, Y, Z ** 10, log=True, colorbar=True)
    viz.plot_cross_section(X, Y, Z, log=False)


def test_plot_3d_current_loops():

    x = np.cos(np.linspace(0, 2 * np.pi))
    y = np.sin(np.linspace(0, 2 * np.pi))
    z = np.zeros_like(y)
    loops = [np.array([x, y, z]).T, 2 * np.array([x, y, z]).T]

    viz.plot_3d_current_loops(loops, tube_radius=0.01)
    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors=None)
