"""
One-element fields
======================

"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh, vector_glyphs

pv.set_plot_theme("document")


import numpy as np

from bfieldtools.mesh_magnetics import (
    scalar_potential_coupling,
    vector_potential_coupling,
)
from bfieldtools.mesh_magnetics import (
    magnetic_field_coupling,
    magnetic_field_coupling_analytic,
)

import trimesh

# Define the element
x = np.sin(np.pi / 6)
y = np.cos(np.pi / 6)
points0 = np.array(
    [[0, 0, 0], [1, 0, 0], [x, y, 0], [-x, y, 0], [-1, 0, 0], [-x, -y, 0], [x, -y, 0]]
)

tris0 = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 6], [0, 6, 1]])
mesh = trimesh.Trimesh(points0, tris0)
scalars = np.zeros(7)
scalars[0] = 1

#%%
#%% Plot element
def plot_element(plotter):
    s1 = plotter.add_mesh(
        triangle_mesh(*points0.T, tris0),
        scalars=scalars,
        cmap="viridis",
        show_scalar_bar=False,
    )
    s2 = plotter.add_mesh(
        triangle_mesh(*points0.T, tris0),
        scalars=scalars,
        cmap="viridis",
        show_scalar_bar=False,
    )
    mesh_actor(s2).mapper.dataset = mesh_actor(s2).mapper.dataset.contour(10)
    mesh_actor(s2).mapper.scalar_range = np.array([0.0, 1.0])
    mesh_actor(s2).mapper.scalar_visibility = False
    mesh_actor(s2).prop.render_lines_as_tubes = True
    mesh_actor(s2).prop.line_width = 3.0


plotter = pv.Plotter()
plot_element(plotter)

#%%
#%% Scalar potential
points = np.array([[0.01, 1, 1], [0.01, 1, -1], [0.01, -1, -1], [0.01, -1, 1]]) * 2
tris = np.array([[0, 1, 2], [2, 3, 0]])
mesh2 = trimesh.Trimesh(points, tris)
for ii in range(7):
    mesh2 = mesh2.subdivide()

U = scalar_potential_coupling(mesh, mesh2.vertices, multiply_coeff=True) @ scalars
plotter = pv.Plotter()
s3 = plotter.add_mesh(
    triangle_mesh(*mesh2.vertices.T, mesh2.faces),
    scalars=U,
    cmap="bwr",
    show_scalar_bar=False,
)
mesh_actor(s3).mapper.dataset = mesh_actor(s3).mapper.dataset.contour(10)
mesh_actor(s3).prop.render_lines_as_tubes = True
mesh_actor(s3).prop.line_width = 3.0
plotter.view_yz()
plot_element(plotter)

#%%
#%% Vector potential
points = np.array([[1, 1, 0.01], [1, -1, 0.01], [-1, -1, 0.01], [-1, 1, 0.01]]) * 2
tris = np.array([[0, 1, 2], [2, 3, 0]])
mesh3 = trimesh.Trimesh(points, tris)
for ii in range(5):
    mesh3 = mesh3.subdivide()
A = vector_potential_coupling(mesh, mesh3.vertices) @ scalars
plotter = pv.Plotter()
vectors = plotter.add_mesh(
    vector_glyphs(*mesh3.vertices.T, *A.T, mode="2ddash"),
    color=(0, 0, 1),
    show_scalar_bar=False,
)
mesh_actor(vectors).prop.render_lines_as_tubes = True
mesh_actor(vectors).prop.line_width = 3.0
plot_element(plotter)

#%%
#%% Magnetic field and its magnitude
from bfieldtools.viz import plot_data_on_vertices

points = (
    np.array([[0.0, 1, 1.001], [0.0, 1, -1], [0.0, -1, -1], [0.0, -1, 1.001]]) * 1.1
)
tris = np.array([[0, 1, 2], [2, 3, 0]])
mesh2 = trimesh.Trimesh(points, tris)
for ii in range(6):
    mesh2 = mesh2.subdivide()

# B0 = magnetic_field_coupling(mesh, mesh2.vertices) @ scalars
B1 = magnetic_field_coupling_analytic(mesh, mesh2.vertices) @ scalars
plotter = plot_data_on_vertices(
    mesh2, np.linalg.norm(B1, axis=1), ncolors=32, colormap="viridis", vmax=1.5e-6
)
vectors = plotter.add_mesh(
    vector_glyphs(*mesh2.vertices.T, *B1.T, scale_factor=50000.0, mode="arrow"),
    color=(1, 0, 1),
    show_scalar_bar=False,
)
plot_element(plotter)
plotter.view_yz()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
