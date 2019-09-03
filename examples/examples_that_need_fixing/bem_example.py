#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 13:37:26 2019

@author: makinea1
"""

import numpy as np
import trimesh
from mayavi import mlab

import pkg_resources
#import sys
#path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
#if path not in sys.path:
#    sys.path.insert(0,path)

from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.integrals import triangle_potential_dipole_linear
from bfieldtools.utils import assemble_matrix
from bfieldtools.magnetic_field_mesh import compute_U

# Load cube representing perfect magnetic shield
file_obj = file_obj=pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/cube.stl')
cube = trimesh.load(file_obj, process=True)
# Center the cube
cube.vertices -= cube.vertices.mean(axis=0)

#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/10x10_plane.obj')
coilmesh = trimesh.load(file_obj, process=False)
# Shrink a little bit and wrap for some additional computations
coilmesh.vertices *= 0.5
coil = MeshWrapper(mesh_obj = coilmesh)

mlab.triangular_mesh(*cube.vertices.T, cube.faces, color=(1,1,1),
                     opacity=0.5)

# Simple constant stream function
weights = np.zeros(coilmesh.vertices.shape[0])
weights[coil.inner_verts] = 1
mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=weights)

############################################
#%% Calculate primary potential matrix
P_prim = compute_U(coilmesh, cube.vertices)

# Plot the resulting primary potential
mlab.figure()
mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=P_prim @ weights,
                     opacity=1.0)

##################################################
#%% Calculate linear collocation BEM matrix
P_bem = compute_U(cube, cube.vertices)

# Recalculate diag elements according to de Munck paper
# Matrix misses one rank
for diag_index in range(P_bem.shape[0]):
    P_bem[diag_index, diag_index] = 0
    P_bem[diag_index, diag_index] = -P_bem[diag_index, :].sum()

# Make it invertible by rank-one update (sets potential of constant dipole layer)
P_bem += np.ones(P_bem.shape)/P_bem.shape[0]

####################################################################
#%% Solve equivalent stream function for the perfect linear mu-metal layer
weights_2 =  np.linalg.solve(P_bem, P_prim @ weights)

#%% Plot the result
mlab.figure()
s1 = mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=weights_2,
                     opacity=1.0)
s1.enable_contours = True
s1.contour.number_of_contours = 30
s2 = mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=weights)
s2.enable_contours = True
s2.contour.number_of_contours = 2

#%%
from bfieldtools.magnetic_field_mesh import compute_C
points = cube.vertices
fpoints1 = points + cube.vertex_normals*0.01

C1_cube = np.moveaxis(compute_C(cube, fpoints1),1,2)
C1_coil = np.moveaxis(compute_C(coilmesh, fpoints1),1,2)

#%% Test field
mlab.figure()
mlab.quiver3d(*points.T, *(-C1_cube @ weights_2 + C1_coil @ weights).T, color=(1,0,0), mode='arrow')
mlab.quiver3d(*points.T, *(C1_cube @ weights_2 + C1_coil @ weights).T, color=(0,0,1), mode='arrow')



