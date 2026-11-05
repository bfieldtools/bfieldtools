"""
Gradient computation example
============================

Minimal example visualizing the gradient of a scalar function on a mesh

"""

import pyvista as pv
from bfieldtools._example_plots import vector_glyphs

pv.set_plot_theme("document")


import numpy as np
import trimesh

from bfieldtools.viz import plot_data_on_vertices
from bfieldtools.mesh_calculus import gradient
from bfieldtools import resources

# Load simple plane mesh that is centered on the origin
file_obj = resources.files("bfieldtools", "example_meshes/10x10_plane.obj")
planemesh = trimesh.load(file_obj, process=False)

# Generate a simple scalar function
r = np.linalg.norm(planemesh.vertices, axis=1)
vals = np.exp(-0.5 * (r / r.max()))

# triangle centers for plotting
tri_centers = planemesh.vertices[planemesh.faces].mean(axis=1).T
tri_centers[1] += 0.1

#%%
#  Calculate the gradient (e.g., flow from potential)
g = gradient(vals, planemesh, rotated=False)

# Plot function and its gradient as arrows

plotter = scene = pv.Plotter(window_size=(800, 800))

plotter = plot_data_on_vertices(planemesh, vals, ncolors=15, figure=scene)
vecs = plotter.add_mesh(
    vector_glyphs(*tri_centers, *g, scale_factor=5, mode="arrow"),
    color=(1, 1, 1),
    show_scalar_bar=False,
)

#%%
# The same but rotated (e.g. current density from a stream function)
g = gradient(vals, planemesh, rotated=True)

# Plot function and its gradient as arrows
plotter = scene = pv.Plotter(window_size=(800, 800))

plotter = plot_data_on_vertices(planemesh, vals, ncolors=15, figure=scene)
vecs = plotter.add_mesh(
    vector_glyphs(*tri_centers, *g, scale_factor=5, mode="arrow"),
    color=(1, 1, 1),
    show_scalar_bar=False,
)

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
