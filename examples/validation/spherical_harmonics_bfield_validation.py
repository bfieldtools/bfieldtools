"""
Spherical harmonics B-field computation validation
==================================================



"""

import numpy as np
import trimesh
from mayavi import mlab

from bfieldtools.mesh_magnetics import magnetic_field_coupling
from bfieldtools.mesh_conductor import MeshConductor

from bfieldtools.sphtools import compute_sphcoeffs_mesh
from bfieldtools import sphtools


import pkg_resources

# Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename(
    "bfieldtools", "example_meshes/10x10_plane.obj"
)
coilmesh = trimesh.load(file_obj, process=False)
coil = MeshConductor(mesh_obj=coilmesh)

coil.mesh.vertices += np.array([0, -1, 0])
weights = np.zeros(coilmesh.vertices.shape[0])
weights[coil.inner_vertices] = 1

test_points = coilmesh.vertices.copy()
test_points[:, 1] = 0

lmax = 9


sph_C = compute_sphcoeffs_mesh(coil.mesh, lmax)

alms = sph_C[0] @ weights
blms = sph_C[1] @ weights

alms = np.zeros_like(alms)

B0 = (magnetic_field_coupling(coilmesh, test_points) @ weights).T
B1 = sphtools.field(test_points, alms, blms, lmax).T

mlab.figure(bgcolor=(1, 1, 1))

s = mlab.triangular_mesh(
    *coilmesh.vertices.T, coilmesh.faces, scalars=weights, colormap="viridis"
)
s.enable_contours = True
s.actor.property.render_lines_as_tubes = True
s.actor.property.line_width = 3.0

mlab.quiver3d(
    *test_points.T, *B0, color=(1, 0, 0), scale_factor=0.5e7, vmin=0, vmax=2e-7
)
mlab.quiver3d(
    *test_points.T, *B1, color=(0, 0, 1), scale_factor=0.5e7, vmin=0, vmax=2e-7
)
s.scene.isometric_view()


print(
    "Relative RMS error", np.sqrt(np.mean((B1 - B0) ** 2)) / np.sqrt(np.mean((B0) ** 2))
)
