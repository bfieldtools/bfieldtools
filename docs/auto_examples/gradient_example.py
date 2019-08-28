'''
Gradient computation example
============================

Minimal example visualizing the gradient of a scalar function on a mesh

'''

import numpy as np
from mayavi import mlab
import trimesh

from bfieldtools.laplacian_mesh import gradient
import pkg_resources

#Load simple plane mesh that is centered on the origin
file_obj = file_obj=pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/10x10_plane_hires.obj')
planemesh = trimesh.load(file_obj, process=False)

# Generate a simple scalar function
r = np.linalg.norm(planemesh.vertices, axis=1)
vals = np.exp(-0.5*(r/r.max()))

# triangle centers for plotting
tri_centers = planemesh.vertices[planemesh.faces].mean(axis=1).T

###############################################################################
#  Calculate the gradient (e.g., flow from potential)
g = gradient(vals, planemesh, rotated=False)

# Plot function and its gradient as arrows
mlab.triangular_mesh(*planemesh.vertices.T, planemesh.faces, scalars=vals)
mlab.quiver3d(*tri_centers, *g, colormap='viridis')

###############################################################################
# The same but rotated (e.g. current density from a stream function)
g = gradient(vals, planemesh, rotated=True)

# Plot function and its gradient as arrows
mlab.figure()
mlab.triangular_mesh(*planemesh.vertices.T, planemesh.faces, scalars=vals)
mlab.quiver3d(*tri_centers, *g, colormap='viridis')
