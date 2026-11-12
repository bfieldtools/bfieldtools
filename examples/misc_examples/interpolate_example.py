"""
Interpolate stream function
===========================

Minimal example showing how to subdivide a mesh and interpolate a scalar function
defined on that mesh to match.

"""

import pyvista as pv
from bfieldtools._example_plots import triangle_mesh, point_cloud

pv.set_plot_theme("document")


import trimesh

import numpy as np
from scipy.linalg import eigh

from bfieldtools.mesh_calculus import laplacian_matrix, mass_matrix
from bfieldtools import utils

from bfieldtools import resources


#%%
# Load a simple mesh and compute an example scalar function on it.
# In this case, the scalar function is an eigenvector of a generalized eigenvalue decomposition

mesh = trimesh.load(resources.files("bfieldtools", "example_meshes/10x10_plane.obj"))

boundaries, inner_verts = utils.find_mesh_boundaries(mesh)

L = laplacian_matrix(mesh)
M = mass_matrix(mesh)

u, v = eigh(
    -L.todense()[inner_verts][:, inner_verts], M.todense()[inner_verts][:, inner_verts]
)

scalars = np.zeros(mesh.vertices.shape[0])
scalars[inner_verts] = v[:, 12]


original_scalars = scalars.copy()
original_mesh = mesh.copy()

#%%
# Plot original scalars and mesh

plotter = scene = pv.Plotter(window_size=(800, 800))

plotter.add_mesh(
    triangle_mesh(*original_mesh.vertices.T, original_mesh.faces),
    scalars=original_scalars,
    style="wireframe",
    show_scalar_bar=False,
)


#%%
# Now, interpolate scalars


ug = triangle_mesh(*mesh.vertices.T, mesh.faces)
ug["scalars"] = scalars


mesh = original_mesh.subdivide().subdivide()
scalars = point_cloud(*mesh.vertices.T).sample(ug)["scalars"]


#%%
# Plot subdivided mesh and interpolated scalars

plotter = scene = pv.Plotter(window_size=(800, 800))

plotter.add_mesh(
    triangle_mesh(*mesh.vertices.T, mesh.faces),
    scalars=scalars,
    style="wireframe",
    show_scalar_bar=False,
)

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
