"""
Analytic B-field computation
==================================================
Validation of analytic mesh operator for magnetic field computation.
"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh, vector_glyphs, point_cloud

pv.set_plot_theme("document")


import numpy as np
import trimesh
import matplotlib.pyplot as plt

from bfieldtools.mesh_calculus import gradient
from bfieldtools.mesh_magnetics import (
    magnetic_field_coupling,
    magnetic_field_coupling_analytic,
)
from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools import resources


# Load simple plane mesh that is centered on the origin
file_obj = resources.files("bfieldtools", "example_meshes/10x10_plane.obj")
coilmesh = trimesh.load(file_obj, process=False)
coil = MeshConductor(mesh_obj=coilmesh)
weights = np.zeros(coilmesh.vertices.shape[0])
weights[coil.inner_vertices] = 1

test_points = coilmesh.vertices + np.array([0, 1, 0])

B0 = magnetic_field_coupling(coilmesh, test_points) @ weights
B1 = magnetic_field_coupling_analytic(coilmesh, test_points) @ weights


plotter = f = pv.Plotter(window_size=(800, 800))

s = plotter.add_mesh(
    triangle_mesh(*coilmesh.vertices.T, coilmesh.faces),
    scalars=weights,
    cmap="viridis",
    show_scalar_bar=False,
)
mesh_actor(s).mapper.dataset = mesh_actor(s).mapper.dataset.contour(10)
mesh_actor(s).prop.render_lines_as_tubes = True
mesh_actor(s).prop.line_width = 3.0

plotter.add_mesh(
    vector_glyphs(*test_points.T, *B0.T), color=(1, 0, 0), show_scalar_bar=False
)
plotter.add_mesh(
    vector_glyphs(*test_points.T, *B1.T), color=(0, 0, 1), show_scalar_bar=False
)

print(
    "Relative RMS error", np.sqrt(np.mean((B1 - B0) ** 2)) / np.sqrt(np.mean((B0) ** 2))
)

#%% Test against analytic formula
# Load simple plane mesh that is centered on the origin
file_obj = resources.files("bfieldtools", "example_meshes/unit_disc.stl")
discmesh = trimesh.load(file_obj, process=True)
for ii in range(3):
    discmesh = discmesh.subdivide()
disc = MeshConductor(mesh_obj=discmesh)
weights = np.zeros(discmesh.vertices.shape[0])
weights[disc.inner_vertices] = 1
plotter = pv.Plotter()
s = plotter.add_mesh(
    triangle_mesh(*discmesh.vertices.T, discmesh.faces),
    scalars=weights,
    cmap="viridis",
    show_scalar_bar=False,
)
g = gradient(weights, discmesh, rotated=True)
plotter.add_mesh(
    vector_glyphs(*discmesh.vertices[discmesh.faces].mean(axis=1).T, *g),
    show_scalar_bar=False,
)

test_points = np.zeros((100, 3))
test_points[:, 2] = np.linspace(0.0, 5, 100)
plotter.add_mesh(
    point_cloud(*test_points.T).glyph(scale=False, geom=pv.Sphere(radius=0.1 / 2)),
    show_scalar_bar=False,
)

# Bfield for 1 Ampere current
B0 = magnetic_field_coupling(discmesh, test_points) @ weights
B1 = magnetic_field_coupling_analytic(discmesh, test_points) @ weights

# Analytic formula for unit disc
plt.plot(1e-7 * 2 * np.pi / (np.sqrt(test_points[:, 2] ** 2 + 1) ** 3))
# Field from the mesh
plt.plot(np.linalg.norm(B0, axis=1))
plt.plot(np.linalg.norm(B1, axis=1))

plt.legend(("Analytic", "Quadrature mesh", "Analytic mesh"))
plt.xlabel("Distance, z [m]")
plt.ylabel("B [T]")

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
