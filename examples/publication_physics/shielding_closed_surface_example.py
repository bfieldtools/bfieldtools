#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 09:25:44 2019

@author: makinea1

Calculate and plot field of a closed shielded current
"""

# import sys
# path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
# if path in sys.path:
#    sys.path.insert(0, path)


import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

from bfieldtools.conductor import Conductor
from bfieldtools.mesh_magnetics import magnetic_field_coupling as compute_C
from bfieldtools.mesh_magnetics import (
    magnetic_field_coupling_analytic as compute_C_analytic,
)
from bfieldtools.mesh_magnetics import scalar_potential_coupling as compute_U
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.mesh_properties import mutual_inductance_matrix
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

from bfieldtools.sphtools import compute_sphcoeffs_mesh

from trimesh.creation import icosphere


mesh1 = icosphere(3, 0.5)
mesh2 = icosphere(3, 1.0)

current1 = Conductor(mesh_obj=mesh1)
current2 = Conductor(mesh_obj=mesh2)

M22 = current2.inductance
M22 += np.ones_like(M22) / M22.shape[0]  # zero mean
M21 = mutual_inductance_matrix(mesh2, mesh1)

x = y = np.linspace(-1.01, 1.01, 150)
X, Y = np.meshgrid(x, y, indexing="ij")
points = np.zeros((X.flatten().shape[0], 3))
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

CB10 = compute_C_analytic(mesh1, points)
CB20 = compute_C_analytic(mesh2, points)

CU10 = compute_U(mesh1, points)
CU20 = compute_U(mesh2, points)

#%%
I1 = mesh1.vertices[:, 0]
# I1 = mesh1.vertices[:,0]**2 - mesh1.vertices[:,1]**2  #- mesh1.vertices[:,2]**2
I2 = -np.linalg.solve(M22, M21 @ I1)

s = mlab.triangular_mesh(*mesh1.vertices.T, mesh1.faces, scalars=I1)
s.enable_contours = True
s = mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=I2)
s.enable_contours = True

B1 = CB10 @ I1
B2 = CB20 @ I2

U1 = CU10 @ I1
U2 = CU20 @ I2
#%% Plot
B = (B1.T + B2.T)[:2].reshape(2, x.shape[0], y.shape[0])
lw = np.sqrt(B[0] ** 2 + B[1] ** 2)
lw = 2 * lw / np.max(lw)
xx = np.linspace(-1, 1, 16)
seed_points = 0.51 * np.array([xx, -np.sqrt(1 - xx ** 2)])
seed_points = np.hstack([seed_points, (0.51 * np.array([xx, np.sqrt(1 - xx ** 2)]))])
# plt.streamplot(x,y, B[1], B[0], density=2, linewidth=lw, color='k',
#               start_points=seed_points.T, integration_direction='both')

#%%
U = (U1 + U2).reshape(x.shape[0], y.shape[0])
U /= np.max(U)
plt.figure()
plt.imshow(
    U,
    vmin=-1.0,
    vmax=1.0,
    cmap="seismic",
    interpolation="bicubic",
    extent=(x.min(), x.max(), y.min(), y.max()),
)
plt.streamplot(
    x,
    y,
    B[1],
    B[0],
    density=2,
    linewidth=lw,
    color="k",
    start_points=seed_points.T,
    integration_direction="both",
)

aa = np.cos(np.linspace(0, 2 * np.pi, 100))
bb = np.sin(np.linspace(0, 2 * np.pi, 100))
plt.plot(0.5 * aa, 0.5 * bb, linewidth=3, color="gray")
plt.plot(1.0 * aa, 1.0 * bb, linewidth=3, color="gray")

plt.xticks([])
plt.yticks([])
