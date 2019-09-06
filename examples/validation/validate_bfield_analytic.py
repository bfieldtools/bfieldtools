#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 11:26:17 2019

@author: makinea1
"""

import numpy as np
import trimesh
from timeit import timeit
from mayavi import mlab
import matplotlib.pyplot as plt

import pkg_resources
import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0,path)

from bfieldtools.laplacian_mesh import gradient
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.magnetic_field_mesh import compute_C_analytic
from bfieldtools.mesh_class import MeshWrapper


#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/10x10_plane.obj')
coilmesh = trimesh.load(file_obj, process=False)
coil = MeshWrapper(mesh_obj = coilmesh)
weights = np.zeros(coilmesh.vertices.shape[0])
weights[coil.inner_verts] = 1

test_points = coilmesh.vertices + np.array([0,1,0])

B0 = np.moveaxis(compute_C(coilmesh, test_points), 2, 0) @ weights
B1 = compute_C_analytic(coilmesh, test_points) @ weights


s = mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces,
                         scalars=weights, colormap='viridis')
s.enable_contours = True
s.actor.property.render_lines_as_tubes = True
s.actor.property.line_width = 3.0

mlab.quiver3d(*test_points.T, *B0, color=(1,0,0))
mlab.quiver3d(*test_points.T, *B1, color=(0,0,1))

np.sum((B1-B0)**2)

#%% Test against

#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/unit_disc.stl')
discmesh = trimesh.load(file_obj, process=True)
for ii in range(3):
    discmesh = discmesh.subdivide()
disc = MeshWrapper(mesh_obj = discmesh)
weights = np.zeros(discmesh.vertices.shape[0])
weights[disc.inner_verts] = 1
mlab.figure()
s = mlab.triangular_mesh(*discmesh.vertices.T, discmesh.faces,
                         scalars=weights, colormap='viridis')
g = gradient(weights, discmesh, rotated=True)
mlab.quiver3d(*discmesh.vertices[discmesh.faces].mean(axis=1).T, *g)

test_points = np.zeros((100, 3))
test_points[:, 2] = np.linspace(0.0, 5, 100)
mlab.points3d(*test_points.T, scale_factor=0.1)

B0 = np.moveaxis(compute_C(discmesh, test_points), 2, 0) @ weights
B1 = compute_C_analytic(discmesh, test_points) @ weights

plt.plot(1e-7*2*np.pi/(np.sqrt(test_points[:,2]**2 + 1)**3))
plt.plot(np.linalg.norm(B0, axis=0))
plt.plot(np.linalg.norm(B1, axis=0))

plt.legend(('Analytic', 'Analytic mesh', 'Quadrature mesh'))


