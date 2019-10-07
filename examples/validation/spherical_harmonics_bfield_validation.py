'''
Spherical harmonics B-field computation validation
==================================================



'''

import numpy as np
import trimesh
from timeit import timeit
from mayavi import mlab
import matplotlib.pyplot as plt

from bfieldtools.laplacian_mesh import gradient
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.magnetic_field_mesh import compute_C_analytic
from bfieldtools.mesh_class import MeshWrapper

from bfieldtools.sphtools import compute_sphcoeffs_mesh, sphbasis


import pkg_resources

#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/10x10_plane.obj')
coilmesh = trimesh.load(file_obj, process=False)
coil = MeshWrapper(mesh_obj = coilmesh)

coil.mesh.vertices += np.array([0,-1,0])
weights = np.zeros(coilmesh.vertices.shape[0])
weights[coil.inner_verts] = 1

test_points = coilmesh.vertices + np.array([0,1,0])

lmax = 10

sph = sphbasis(20)

sph_C = compute_sphcoeffs_mesh(coil.mesh, lmax)

alms = sph_C[0] @ weights
blms = sph_C[1] @ weights

blms = np.zeros_like(alms)

B0 = np.moveaxis(compute_C(coilmesh, test_points), 2, 0) @ weights
B1 = sph.field(test_points, alms, blms, lmax).T



s = mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces,
                         scalars=weights, colormap='viridis')
s.enable_contours = True
s.actor.property.render_lines_as_tubes = True
s.actor.property.line_width = 3.0

mlab.quiver3d(*test_points.T, *B0, color=(1,0,0))
mlab.quiver3d(*test_points.T, *B1, color=(0,0,1))

print('Relative RMS error',  np.sqrt(np.mean((B1-B0)**2))/np.sqrt(np.mean((B0)**2)))