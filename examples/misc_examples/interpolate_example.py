"""
Interpolate stream function
===========================

Minimal example showing how to subdivide a mesh and interpolate a scalar function
defined on that mesh to match.

"""

from tvtk.api import tvtk
from mayavi import mlab
import trimesh

import numpy as np
from scipy.linalg import eigh

from bfieldtools.mesh_calculus import laplacian_matrix, mass_matrix
from bfieldtools import utils

import pkg_resources


#%%
# Load a simple mesh and compute an example scalar function on it.
# In this case, the scalar function is an eigenvector of a generalized eigenvalue decomposition

mesh = trimesh.load(
    pkg_resources.resource_filename("bfieldtools", "example_meshes/10x10_plane.obj")
)

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

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))

mlab.triangular_mesh(
    *original_mesh.vertices.T,
    original_mesh.faces,
    scalars=original_scalars,
    representation="wireframe"
)


#%%
# Now, interpolate scalars


ug = tvtk.UnstructuredGrid(points=mesh.vertices)

ug.set_cells(tvtk.Triangle().cell_type, mesh.faces)
ug.point_data.scalars = scalars
ug.point_data.scalars.name = "scalars"


mesh = original_mesh.subdivide().subdivide()
scalars = mlab.pipeline.probe_data(ug, *mesh.vertices.T)


#%%
# Plot subdivided mesh and interpolated scalars

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))

mlab.triangular_mesh(
    *mesh.vertices.T, mesh.faces, scalars=scalars, representation="wireframe"
)
