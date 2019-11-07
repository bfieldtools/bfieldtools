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
from bfieldtools.suhtools import suhbasis

# Load cube representing perfect magnetic shield
file_obj = file_obj=pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/closed_cylinder.stl')
cube = trimesh.load(file_obj, process=True)
# Center the cube
cube.vertices -= cube.vertices.mean(axis=0)
cube = cube.subdivide()

#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/cylinder.stl')
coilmesh = trimesh.load(file_obj, process=True)
coilmesh=coilmesh.subdivide()
# Shrink a little bit and wrap for some additional computations
coilmesh.vertices = np.roll(coilmesh.vertices, 1, 1)
coilmesh.vertices *= 0.22
coilmesh.vertices[:, 0] *= 0.5
coilmesh.vertices[:, 0] += 0.2
coilmesh.vertices[:, 2] += 0.05
coil = MeshWrapper(mesh_obj = coilmesh)
1
#%%
mlab.figure(bgcolor=(1,1,1))
mlab.triangular_mesh(*cube.vertices.T, cube.faces, color=(1,1,1),
                     opacity=0.5)

sb= suhbasis(coilmesh, Nc=5, closed_mesh=False)
weights = sb.basis[:,4]
# Simple constant stream function
#weights = np.zeros(coilmesh.vertices.shape[0])
#weights[coil.inner_verts] = 1
#weights = coilmesh.vertices[:, 0]**2 - coilmesh.vertices[:, 2]**2

s2 = mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=weights)
s2.enable_contours = True
s2.contour.number_of_contours = 10
############################################
#%% Calculate primary potential matrix
P_prim = compute_U(coilmesh, cube.vertices)

#%% Plot the resulting primary potential
mlab.figure()
s =mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=P_prim @ weights,
                     opacity=1.0)
s.enable_contours = True
s.contour.filled_contours = True
s.contour.number_of_contours = 30

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
#P_bem += 2*np.pi*np.eye(P_bem.shape[0])
####################################################################
#%% Solve equivalent stream function for the perfect linear mu-metal layer
weights_2 =  np.linalg.solve(-P_bem, P_prim @ weights)

#%% Plot the result
mlab.figure()
s1 = mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=weights_2,
                     opacity=1.0)
s1.enable_contours = True
#s1.contour.filled_contours = True
s1.contour.number_of_contours = 30
#s2 = mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=weights)
#s2.enable_contours = True
#s2.contour.number_of_contours = 2

#%%
from bfieldtools.magnetic_field_mesh import compute_C_analytic
#points = np.mean(cube.vertices[cube.faces], axis=1)
#fpoints1 = points #+ cube.face_normals*0.01

#C1_cube = compute_C_analytic(cube, fpoints1)
#C1_coil = compute_C_analytic(coilmesh, fpoints1)

#%% Test field
#mlab.figure()
#mlab.quiver3d(*points.T, *(-C1_cube @ weights_2 + C1_coil @ weights).T, color=(1,0,0), mode='arrow')
#mlab.quiver3d(*points.T, *(C1_cube @ weights_2 + C1_coil @ weights).T, color=(0,0,1), mode='arrow')

#%%
#points = np.mean(coilmesh.vertices[coilmesh.faces], axis=1)
N= 30
points = np.zeros((N*N, 3))
Y,Z = np.meshgrid(np.linspace(-0.1, 0.1, N), np.linspace(-0.1, 0.1, N), indexing='ij')
points[:, 1] = Y.flatten()
points[:, 2] = Z.flatten()
points[:, 0] += coilmesh.vertices.mean(axis=0)
points = points[np.linalg.norm(points, axis=1) < 4]
fpoints2 = points #+ cube.face_normals*0.01

C2_cube = compute_C_analytic(cube, fpoints2)
C2_coil = compute_C_analytic(coilmesh, fpoints2)

#%% Test field
mlab.figure()
mlab.quiver3d(*points.T, *(0*C2_cube @ weights_2 + C2_coil @ weights).T, color=(1,0,0), mode='arrow')
mlab.quiver3d(*points.T, *(C2_cube @ weights_2 + C2_coil @ weights).T, color=(0,0,1), mode='arrow')



