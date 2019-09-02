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


# Load cube representing perfect magnetic shield
file_obj = file_obj=pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/cube.stl')
cube = trimesh.load(file_obj, process=False)
# Center the cube
cube.vertices -= cube.vertices.mean(axis=0)

#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/10x10_plane_hires.obj')
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

# Calculate difference vectors
R1 = coilmesh.vertices[coilmesh.faces]
R2 = cube.vertices
RR = R2[:, None, None, :] - R1[None, :, :, :]

############################################
#%% Calculate primary potential matrix in chunks
RRchunks = np.array_split(RR, 10, axis=0)
i0=0;
Pf = np.zeros((R2.shape[0], coilmesh.faces.shape[0], 3))
print('Computing coupling matrix in chunks')
for RRchunk in RRchunks:
    i1 = i0+RRchunk.shape[0]
    Pi = triangle_potential_dipole_linear(RRchunk, coilmesh.face_normals,
                                         coilmesh.area_faces)
    Pf[i0:i1] = Pi
    print(i1/R2.shape[0]*100, '% computed')
    i0=i1
# Accumulate the elements
Pv = np.zeros((R2.shape[0], coilmesh.vertices.shape[0]))
Pv[:, coilmesh.faces[:, 0]] += Pf[:, :, 0]
Pv[:, coilmesh.faces[:, 1]] += Pf[:, :, 1]
Pv[:, coilmesh.faces[:, 2]] += Pf[:, :, 2]

P_prim=Pv.copy()

# Plot the resulting primary potential
mlab.figure()
mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=Pv @ weights,
                     opacity=1.0)

##################################################
#%% Calculate linear collocation BEM matrix in chunks

# Source and eval locations
R1 = cube.vertices[cube.faces]
R2 = cube.vertices

R2chunks = np.array_split(R2, 20, axis=0)
i0=0;
Pf = np.zeros((R2.shape[0], cube.faces.shape[0], 3))
print('Computing BEM matrix in chunks')
for R2chunk in R2chunks:
    RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
    i1 = i0+RRchunk.shape[0]
    Pi = triangle_potential_dipole_linear(RRchunk, cube.face_normals,
                                         cube.area_faces)
    Pf[i0:i1] = Pi
    print((100*i1)//R2.shape[0], '% computed')
    i0=i1

# Accumulate the elements
Pv = np.zeros((R2.shape[0], cube.vertices.shape[0]))
Pv[:, cube.faces[:, 0]] += Pf[:, :, 0]
Pv[:, cube.faces[:, 1]] += Pf[:, :, 1]
Pv[:, cube.faces[:, 2]] += Pf[:, :, 2]

P_bem = Pv

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

#PLot
s1 = mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=weights_2,
                     opacity=1.0)
s1.enable_contours = True
s1.contour.number_of_contours = 30
s2 = mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=weights)
s2.enable_contours = True
s2.contour.number_of_contours = 2

