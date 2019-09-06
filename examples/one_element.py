#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 15:46:47 2019

@author: makinea1
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0,path)

from bfieldtools.integrals import triangle_potential_dipole_linear
from bfieldtools.integrals import omega
from bfieldtools.utils import tri_normals_and_areas
from bfieldtools.laplacian_mesh import gradient
from bfieldtools.magnetic_field_mesh import compute_U, compute_A
from bfieldtools.magnetic_field_mesh import compute_C, compute_C_analytic

import trimesh
from mayavi import mlab

#########################################################
#%% Test potential shape slightly above the surface
x = np.sin(np.pi/6)
y = np.cos(np.pi/6)
points = np.array([[0, 0, 0],
                   [1, 0, 0],
                   [x, y, 0],
                   [-x, y, 0],
                   [-1, 0, 0],
                   [-x, -y, 0],
                   [x, -y, 0]])

tris = np.array([[0,1,2],[0,2,3],[0,3,4],[0,4,5],[0,5,6],[0,6,1]])
mesh = trimesh.Trimesh(points, tris)
scalars = np.zeros(7)
scalars[0] = 1
# Stream function
s1 = mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap='viridis')
# Stream lines
s2 = mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap='viridis')
s2.enable_contours = True
s2.actor.mapper.scalar_range = np.array([0., 1.])
s2.actor.mapper.scalar_visibility = False
s2.actor.property.render_lines_as_tubes = True
s2.actor.property.line_width = 3.0

#%%
points = np.array([[0.01, 1, 1],
                   [0.01, 1, -1],
                   [0.01, -1, -1],
                   [0.01, -1, 1]])*2
tris=np.array([[0,1,2], [2,3,0]])
mesh2 = trimesh.Trimesh(points, tris)
for ii in range(7):
    mesh2 =mesh2.subdivide()

U = compute_U(mesh, mesh2.vertices) @ scalars

s3= mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=U, colormap='bwr')
s3.enable_contours = True
s3.contour.minimum_contour = -5.2e-07
s3.contour.maximum_contour = 5.2e-07
s3.actor.property.render_lines_as_tubes = True
s3.actor.property.line_width = 3.0

#%%
if False:
    points = np.array([[1, 1, -0.01],
                       [1, -1, -0.01],
                       [-1, -1, -0.01],
                       [-1, 1, -0.01]])*2
    tris=np.array([[0,1,2], [2,3,0]])
    mesh3 = trimesh.Trimesh(points, tris)
    for ii in range(5):
        mesh3 =mesh3.subdivide()
    A = compute_A(mesh, mesh3.vertices) @ scalars
    vectors = mlab.quiver3d(*mesh3.vertices.T, *A, mode='2ddash', color=(0,0,1))
    vectors.glyph.glyph_source.glyph_position = 'center'
    vectors.actor.property.render_lines_as_tubes = True
    vectors.actor.property.line_width = 3.0
#%%
#B0 = np.moveaxis(compute_C(mesh, mesh2.vertices), 2, 0) @ scalars
points = np.array([[0.001, 1, 1],
                   [0.001, 1, -1],
                   [0.001, -1, -1],
                   [0.001, -1, 1]])*2 + 0.001
tris=np.array([[0,1,2], [2,3,0]])
mesh2 = trimesh.Trimesh(points, tris)
for ii in range(6):
    mesh2 =mesh2.subdivide()
B1 = compute_C_analytic(mesh, mesh2.vertices) @ scalars
vectors = mlab.quiver3d(*mesh2.vertices.T, *B1, mode='arrow', color=(1,0,0))
vectors.glyph.glyph_source.glyph_position = 'center'
#vectors.actor.property.render_lines_as_tubes = True
#vectors.actor.property.line_width = 3.0


