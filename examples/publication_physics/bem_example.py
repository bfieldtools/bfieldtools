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
import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0,path)

from bfieldtools.mesh_class import Conductor
from bfieldtools.integrals import triangle_potential_dipole_linear
from bfieldtools.utils import assemble_matrix
from bfieldtools.magnetic_field_mesh import compute_U
from bfieldtools.suhtools import suhbasis

# Load cube representing perfect magnetic shield
file_obj = file_obj=pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/closed_cylinder_remeshed.stl')
cube = trimesh.load(file_obj, process=True)
# Center the cube
cube.vertices -= cube.vertices.mean(axis=0)
#cube = cube.subdivide()

#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/cylinder.stl')
coilmesh = trimesh.load(file_obj, process=True)
#coilmesh=coilmesh.subdivide()
# Shrink a little bit and wrap for some additional computations
coilmesh.vertices = np.roll(coilmesh.vertices, 1, 1)
coilmesh.vertices *= 0.2
coilmesh.vertices[:, 0] *= 0.5
coilmesh.vertices[:, 0] += 0.25
coilmesh.vertices[:, 1:2] -= coilmesh.vertices[:, 1:2].mean(axis=0)
coilmesh.vertices[:, 1:2] += cube.vertices[:, 1:2].mean(axis=0)
coil = Conductor(mesh_obj = coilmesh)

#%%
fig = mlab.figure(bgcolor=(1,1,1))
s0 = mlab.triangular_mesh(*cube.vertices.T, cube.faces, color=(0.5,0.5,0.5),
                     opacity=0.3)
s0.actor.property.backface_culling = False

sb= suhbasis(coilmesh, Nc=9, closed_mesh=False)
weights = sb.basis[:,1]
#weights = abs(weights)
# Simple constant stream function
#weights = np.zeros(coilmesh.vertices.shape[0])
#weights[coil.inner_verts] = 1
#weights = coilmesh.vertices[:, 0]**2 - coilmesh.vertices[:, 2]**2

s2 = mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=weights,
                          colormap='RdBu')
s2.enable_contours = True
#s2.contour.filled_contours = True
s2.contour.number_of_contours = 20
s2.actor.property.render_lines_as_tubes = True
s2.actor.property.ambient = 0.2


scene = s0.module_manager
scene.scene.camera.position = [-2.439666317033772, 0.27339154055258547, 4.047442287738252]
scene.scene.camera.focal_point = [0.006026178483206079, -0.004020361968711389, 0.005475058591664084]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.28353553256074365, 0.9530687973708336, 0.10614833608497154]
scene.scene.camera.clipping_range = [2.7484502681121623, 7.241260542714118]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
############################################
#%% Calculate primary potential matrix
# Compute slightly inside
d = 1e-3
P_prim = compute_U(coilmesh, cube.vertices - d*cube.vertex_normals)

#%% Plot the resulting primary potential
mlab.figure()
s =mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=P_prim @ weights,
                     opacity=1.0)
s.enable_contours = True
s.contour.filled_contours = True
s.contour.number_of_contours = 30

##################################################
#%% Calculate linear collocation BEM matrix
P_bem = compute_U(cube, cube.vertices - d*cube.vertex_normals)

####################################################################
#%% Solve equivalent stream function for the perfect linear mu-metal layer
weights_2 =  np.linalg.solve(-P_bem, P_prim @ weights)

#%% Plot the result
fig = mlab.figure(bgcolor=(1,1,1))
s0 = mlab.triangular_mesh(*cube.vertices.T, cube.faces, color=(0.5,0.5,0.5),
                     opacity=0.3)
s0.actor.property.backface_culling = False
s1 = mlab.triangular_mesh(*cube.vertices.T, cube.faces, scalars=weights_2,
                     opacity=1.0, colormap='RdBu',
                     vmin=-max(abs(weights_2)), vmax=max(abs(weights_2)))
s1.enable_contours = True
s1.actor.property.render_lines_as_tubes = True

#s1.contour.filled_contours = True
s1.contour.number_of_contours = 30
s1.actor.property.render_lines_as_tubes = True
s1.actor.property.ambient = 0.2

scene = s1.module_manager
scene.scene.camera.position = [-2.439666317033772, 0.27339154055258547, 4.047442287738252]
scene.scene.camera.focal_point = [0.006026178483206079, -0.004020361968711389, 0.005475058591664084]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.28353553256074365, 0.9530687973708336, 0.10614833608497154]
scene.scene.camera.clipping_range = [2.7484502681121623, 7.241260542714118]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()

#%%
#from bfieldtools.magnetic_field_mesh import compute_C_analytic
#points = np.mean(cube.vertices[cube.faces], axis=1)
#fpoints1 = points - cube.face_normals*d
##points = cube.vertices
##fpoints1 = points - cube.vertex_normals*d
#C1_cube = compute_C_analytic(cube, fpoints1)
#C1_coil = compute_C_analytic(coilmesh, fpoints1)

#%% Test field
#mlab.figure()
#mlab.quiver3d(*points.T, *(-C1_cube @ weights_2 + C1_coil @ weights).T, color=(1,0,0), mode='arrow')
#mlab.quiver3d(*points.T, *(C1_cube @ weights_2 + C1_coil @ weights).T, color=(0,0,1), mode='arrow')

#%%
#N= 30
#points = np.zeros((N*N, 3))
#Y,Z = np.meshgrid(np.linspace(-0.1, 0.1, N), np.linspace(-0.1, 0.1, N), indexing='ij')
#points[:, 0] = Y.flatten()
#points[:, 2] = Z.flatten()
#points += coilmesh.vertices.mean(axis=0)
#points = points[np.linalg.norm(points, axis=1) < 4]
#fpoints2 = points #+ cube.face_normals*0.01
#
#C2_cube = compute_C_analytic(cube, fpoints2)
#C2_coil = compute_C_analytic(coilmesh, fpoints2)
#
##%% Test field
##mlab.figure()
#mlab.quiver3d(*points.T, *(0*C2_cube @ weights_2 + C2_coil @ weights).T, color=(1,0,0), mode='arrow')
#mlab.quiver3d(*points.T, *(C2_cube @ weights_2 + C2_coil @ weights).T, color=(0,0,1), mode='arrow')
#

#%%
N= 80
points = np.zeros((N*N, 3))
w = 0.55
X,Z = np.meshgrid(np.linspace(-w, w, N), np.linspace(-w, w, N), indexing='ij')
points[:, 0] = X.flatten()
points[:, 2] = Z.flatten()
X+= coilmesh.vertices.mean(axis=0)[0]
Z+= coilmesh.vertices.mean(axis=0)[2]
points += coilmesh.vertices.mean(axis=0)
fpoints2 = points
U2_shield = compute_U(cube, fpoints2)
U2_prim = compute_U(coilmesh, fpoints2)

#%%
from bfieldtools.contour import scalar_contour
cc1 = scalar_contour(coilmesh, coilmesh.vertices[:,1], contours= [-0.001])[0]
cc1a = cc1[0]
cc1b = np.vstack(cc1[1:])
cc2 = scalar_contour(cube, cube.vertices[:,1], contours= [-0.001])[0][0]
cc2 = cc2[cc2[:,0] > 0]



import matplotlib.pyplot as plt
u0 = abs(np.sum(U2_shield, axis=1).reshape(N, N)) # Solid angle of the shield, zero outside
u0 /= u0.max()
u0[u0 < 1e-6]  = 0
u1 = (U2_prim @ weights).reshape(N,N)
u2 = (U2_shield @ weights_2).reshape(N,N)*u0
u3= (u1 + u2)*u0

levels = np.linspace(-np.max(abs(u3)), np.max(abs(u3)), 64)
p =plt.contourf(X,Z, u1, levels=levels, cmap='seismic')
plt.plot(cc1a[:,0],cc1a[:,2], linewidth=3, color='gray')
plt.plot(cc1b[:,0],cc1b[:,2], linewidth=3, color='gray')
plt.axis('image')
plt.axis('off')
xlims = p.ax.get_xlim()
plt.figure()
p = plt.contourf(X,Z, u2, levels=levels, cmap='seismic')
plt.plot(cc1a[:,0],cc1a[:,2], linewidth=3, color='gray')
plt.plot(cc1b[:,0],cc1b[:,2], linewidth=3, color='gray')
plt.plot(cc2[:,0],cc2[:,2], linewidth=3, color='gray')
plt.axis('image')
plt.axis('off')
p.ax.set_xlim(xlims)
plt.figure()
p =plt.contourf(X,Z, u3, levels=levels, cmap='seismic')
plt.plot(cc1a[:,0],cc1a[:,2], linewidth=3, color='gray')
plt.plot(cc1b[:,0],cc1b[:,2], linewidth=3, color='gray')
plt.plot(cc2[:,0],cc2[:,2], linewidth=3, color='gray')
plt.axis('image')
plt.axis('off')
p.ax.set_xlim(xlims)
