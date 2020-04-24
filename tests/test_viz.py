# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:09:56 2020

@author: Rasmus Zetter
"""

import pytest
from bfieldtools import viz
from bfieldtools.utils import load_example_mesh

from .test_mesh_conductor import _fake_mesh_conductor, _fake_streamfunction


import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt


def test_mesh_plotting():

    mesh = load_example_mesh("unit_disc")

    v_scalars = np.random.rand(len(mesh.vertices))
    f_scalars = np.random.rand(len(mesh.faces))

    viz.plot_mesh(mesh, figure=None)
    f = mlab.figure()
    viz.plot_mesh(mesh, figure=f)

    viz.plot_data_on_faces(mesh, f_scalars)

    f = mlab.figure()
    viz.plot_data_on_faces(mesh, f_scalars, figure=f, colormap="jet")

    f_scalars = np.random.rand(len(mesh.faces)) - 0.5
    viz.plot_data_on_faces(mesh, f_scalars, figure=f, colorbar=True)

    viz.plot_data_on_vertices(mesh, v_scalars, colorbar=True)
    viz.plot_data_on_vertices(mesh, v_scalars, autoscale=True)

    v_scalars = np.random.rand(len(mesh.vertices)) - 0.5

    viz.plot_data_on_vertices(mesh, v_scalars, colormap="jet")


def test_plot_cross_section():

    X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))
    Z = 3 - (X ** 2 + Y ** 2)

    viz.plot_cross_section(X, Y, Z ** 10, log=True, colorbar=True, cmap=None)
    viz.plot_cross_section(X, Y, Z, log=False, cmap=None)
    viz.plot_cross_section(
        X, Y, Z, log=False, cmap="jet", vmin=-10, vmax=10, contours=False
    )

    X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))
    Z = 3 - 9 * (X ** 2 + Y ** 2)

    viz.plot_cross_section(X, Y, Z, log=False, cmap=None, contours=4)

    fig, axes = plt.subplots(1, 1)
    viz.plot_cross_section(X, Y, Z, log=False, axes=axes, vmin=None, vmax=None)


def test_plot_3d_current_loops():

    x = np.cos(np.linspace(0, 2 * np.pi))
    y = np.sin(np.linspace(0, 2 * np.pi))
    z = np.zeros_like(y)
    loops = [np.array([x, y, z]).T, 2 * np.array([x, y, z]).T]

    viz.plot_3d_current_loops(loops, tube_radius=0.01)
    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors="auto")
    f = mlab.figure()
    viz.plot_3d_current_loops(loops, tube_radius=0.01, figure=f)

    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors=None)
    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors=(1, 0, 0))
    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors=[(1, 0, 0)] * 2)

    # Longest segment is the one closing the loop
    x = np.cos(np.linspace(0, 2 * np.pi))[:-1]
    y = np.sin(np.linspace(0, 2 * np.pi))[:-1]
    z = np.zeros_like(y)
    loops = [np.array([x, y, z]).T, 2 * np.array([x, y, z]).T]
    viz.plot_3d_current_loops(loops, tube_radius=0.01)


@pytest.mark.xfail
def test_plot_3d_current_loops_fail():

    x = np.cos(np.linspace(0, 2 * np.pi))
    y = np.sin(np.linspace(0, 2 * np.pi))
    z = np.zeros_like(y)
    loops = [np.array([x, y, z]).T, 2 * np.array([x, y, z]).T]

    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors="boop")
