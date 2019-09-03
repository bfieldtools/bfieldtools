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
from bfieldtools.magnetic_field_mesh import compute_U

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
s1 = mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap='viridis')

s2 = mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap='viridis')
s2.enable_contours = True
s2.actor.mapper.scalar_range = np.array([0., 1.])
s2.actor.mapper.scalar_visibility = False
s2.actor.property.render_lines_as_tubes = True

points = np.array([[0.01, 1, 1],
                   [0.01, 1, -1],
                   [0.01, -1, -1],
                   [0.01, -1, 1]])*1.5
tris=np.array([[0,1,2], [2,3,0]])
mesh2 = trimesh.Trimesh(points, tris)
for ii in range(7):
    mesh2 =mesh2.subdivide()

U = compute_U(mesh, mesh2.vertices) @ scalars

s3= mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=U, colormap='bwr')
s3.enable_contours = True
s3.actor.mapper.scalar_range = np.array([0., 1.])
s3.actor.property.render_lines_as_tubes = True