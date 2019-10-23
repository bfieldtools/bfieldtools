#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 09:25:44 2019

@author: makinea1

Calculate and plot field of a closed shielded current
"""

#import sys
#path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
#if path in sys.path:
#    sys.path.insert(0, path)


import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh

from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.magnetic_field_mesh import compute_C_analytic
from bfieldtools.magnetic_field_mesh import compute_U
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix_from_A
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops
from bfieldtools.sphtools import compute_sphcoeffs_mesh, sphbasis

from trimesh.creation import box
from trimesh.smoothing import filter_laplacian

import pkg_resources

#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 0.1

#Load simple plane mesh that is centered on the origin
planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

planemesh.apply_scale(scaling_factor)

#Specify coil plane geometry
center_offset = np.array([0, 0, 0]) * scaling_factor
standoff = np.array([0, 4, 0]) * scaling_factor

#Create coil plane pairs
coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                         planemesh.faces, process=False)

coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                     planemesh.faces, process=False)

mesh1 = coil_plus.union(coil_minus)
mesh2 = mesh1.copy()
mesh2.apply_scale(1.4)

current1 = MeshWrapper(mesh_obj = mesh1, fix_normals=False)
current2 = MeshWrapper(mesh_obj = mesh2, fix_normals=False)



#
#mesh1 = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)
#mesh1.vertices *= 0.1
#mesh1.vertices[:,1] += 0.5
#
#mesh2 = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)
#mesh2.vertices *= 0.15
#mesh2.vertices[:,1] += 0.75


current1 = MeshWrapper(mesh_obj=mesh1)
current2 = MeshWrapper(mesh_obj=mesh2)

M11 = current1.inductance
M22 = current2.inductance
# Constrain boundary to zero and consider only inneverts
M11 = M11[current1.inner_verts][:, current1.inner_verts]
M22 = M22[current2.inner_verts][:, current2.inner_verts]
# Add rank-one matrix, so that M22 can be inverted (for zero mean functions)
#M22 += np.ones_like(M22)/M22.shape[0]
#M11 += np.ones_like(M11)/M11.shape[0]
M21 = mutual_inductance_matrix_from_A(mesh2, mesh1)
M21 = M21[current2.inner_verts][:, current1.inner_verts]
# Mapping from I1 to I2, constraining flux through mesh2 to zero
P = -np.linalg.solve(M22, M21)

#A1, temp = compute_sphcoeffs_mesh(mesh1, 7)
#A2, temp = compute_sphcoeffs_mesh(mesh2, 7)

sb = sphbasis(10)
F1 = (sb.basis_fields(mesh1.vertices, 3)[1]*mesh1.vertex_normals).sum(axis=-1)
#F2 = (sb.basis_fields(mesh2.vertices, 3)[0]*mesh2.vertex_normals).sum(axis=-1)

x = y = np.linspace(-0.8, 0.8, 50)
X,Y = np.meshgrid(x, y, indexing='ij')
points = np.zeros((X.flatten().shape[0], 3))
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

CB1 = compute_C_analytic(mesh1, points)
CB2 = compute_C_analytic(mesh2, points)

CU1 = compute_U(mesh1, points)
CU2 = compute_U(mesh2, points)

#%% Specify spherical harmonic and calculate corresponding shielded field
#alpha = np.zeros(A1.shape[0])
#alpha[2] = 1
#alpha[15] = 1
# Minimization of magnetic energy with spherical harmonic constraint
#C = A1 + A2 @ P
M = M11 + M21.T @ P
#G = np.linalg.solve(M, C.T)
#I1 = G @ np.linalg.solve(C @ G, alpha)

#f = F1[2]
f =  np.ones(M11.shape[0])
#f = mesh1.vertex_normals[:,1]
I1inner = np.linalg.solve(M, f)
#I1, res, rr, s = np.linalg.lstsq(C, alpha, rcond=1e-12)
#I1 = mesh1.vertices[:,0]
I2inner = P @ I1inner

I1 = np.zeros(mesh1.vertices.shape[0]); I1[current1.inner_verts] = I1inner
I2 = np.zeros(mesh2.vertices.shape[0]); I2[current2.inner_verts] = I2inner

s = mlab.triangular_mesh(*mesh1.vertices.T, mesh1.faces, scalars=I1)
s.enable_contours=True
s = mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=I2)
s.enable_contours=True

B1 = CB1 @ I1
B2 = CB2 @ I2

U1 = CU1 @ I1
U2 = CU2 @ I2
#%% Plot
B = (B1.T + B2.T)[:2].reshape(2, x.shape[0], y.shape[0])
lw = np.sqrt(B[0]**2 + B[1]**2)
lw = 2*lw/np.max(lw)
xx = np.linspace(-1,1, 16)
seed_points = 0.51*np.array([xx, -np.sqrt(1-xx**2)])
seed_points = np.hstack([seed_points, (0.51*np.array([xx, np.sqrt(1-xx**2)]))])
#plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
#               start_points=seed_points.T, integration_direction='both')
U = (U1 + U2).reshape(x.shape[0], y.shape[0])
U /= np.max(U)
plt.figure()
plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
           extent=(x.min(), x.max(), y.min(), y.max()))
plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
               start_points=seed_points.T, integration_direction='both')

cc1 = scalar_contour(mesh1, mesh1.vertices[:,2], contours= [-0.001])[0][0]
cc2 = scalar_contour(mesh2, mesh2.vertices[:,2], contours= [-0.001])[0][0]

plt.plot(cc1[:,1], cc1[:,0], linewidth=3.0)
plt.plot(cc2[:,1], cc2[:,0], linewidth=3.0)

plt.xticks([])
plt.yticks([])





