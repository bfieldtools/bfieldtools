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
import pyvista as pv
import matplotlib.pyplot as plt


def test_mesh_plotting():

    mesh = load_example_mesh("unit_disc")

    v_scalars = np.random.rand(len(mesh.vertices))
    f_scalars = np.random.rand(len(mesh.faces))

    viz.plot_mesh(mesh, figure=None)
    f = pv.Plotter()
    viz.plot_mesh(mesh, figure=f)

    viz.plot_data_on_faces(mesh, f_scalars)

    f = pv.Plotter()
    viz.plot_data_on_faces(mesh, f_scalars, figure=f, colormap="jet")

    f = pv.Plotter()
    f_scalars = np.random.rand(len(mesh.faces)) - 0.5
    viz.plot_data_on_faces(mesh, f_scalars, figure=f, colorbar=True)

    viz.plot_data_on_vertices(mesh, v_scalars, colorbar=True)
    viz.plot_data_on_vertices(mesh, v_scalars, autoscale=True)

    v_scalars = np.random.rand(len(mesh.vertices)) - 0.5

    viz.plot_data_on_vertices(mesh, v_scalars, colormap="jet")


def test_plot_cross_section():

    X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))
    Z = 3 - (X**2 + Y**2)

    viz.plot_cross_section(X, Y, Z**10, log=True, colorbar=True, cmap=None)
    viz.plot_cross_section(X, Y, Z, log=False, cmap=None)
    viz.plot_cross_section(
        X, Y, Z, log=False, cmap="jet", vmin=-10, vmax=10, contours=False
    )

    X, Y = np.meshgrid(np.linspace(0, 1, 11), np.linspace(0, 1, 11))
    Z = 3 - 9 * (X**2 + Y**2)

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
    f = pv.Plotter()
    viz.plot_3d_current_loops(loops, tube_radius=0.01, figure=f)

    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors=None)
    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors="b")
    viz.plot_3d_current_loops(loops, tube_radius=0.01, colors=["b"] * 2)

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


@pytest.mark.parametrize("on_faces", [False, True])
@pytest.mark.parametrize(
    "vmin,vmax,expected", [(None, 3.0, (-1.0, 3.0)), (-2.0, None, (-2.0, 2.0))]
)
def test_partial_colour_limits(on_faces, vmin, vmax, expected):
    mesh = load_example_mesh("unit_disc")
    plotter = pv.Plotter(off_screen=True)
    try:
        count = len(mesh.faces) if on_faces else len(mesh.vertices)
        data = np.linspace(-1.0, 2.0, count)
        plot = viz.plot_data_on_faces if on_faces else viz.plot_data_on_vertices
        plot(mesh, data, figure=plotter, colormap="bwr", vmin=vmin, vmax=vmax)
        actor = next(a for a in plotter.actors.values() if isinstance(a, pv.Actor))
        np.testing.assert_allclose(actor.mapper.scalar_range, expected)
    finally:
        plotter.close()


def test_current_loops_single_rgb_colour():
    angle = np.linspace(0, 2 * np.pi, 20)
    loop = np.column_stack((np.cos(angle), np.sin(angle), np.zeros_like(angle)))
    plotter = pv.Plotter(off_screen=True)
    try:
        viz.plot_3d_current_loops([loop], figure=plotter, colors=(0.9, 0.8, 0.7))
        actor = next(a for a in plotter.actors.values() if isinstance(a, pv.Actor))
        np.testing.assert_allclose(
            actor.prop.color.float_rgb, (0.9, 0.8, 0.7), atol=1 / 255
        )
    finally:
        plotter.close()
