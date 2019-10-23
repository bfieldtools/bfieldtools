"""
Plot different field patterns the basic integrals

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from mayavi import mlab
import trimesh

path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0,path)

from bfieldtools.integrals import triangle_potential_uniform
from bfieldtools.integrals import triangle_potential_dipole_linear
from bfieldtools.integrals import gamma0
from bfieldtools.integrals import omega
from bfieldtools.utils import tri_normals_and_areas

#########################################################
#%% Create there orthogonal planes
points = np.array([[0,0,0],
                   [1,0.5,0],
                   [0,1,0]]) + 0.15

points[:,2] += 0.1

tris = np.array([[0,1,2]])
p_tris = points[tris]
mesh = trimesh.Trimesh(points, tris, process=False)

# Evaluation points
points2 = np.array([[-0.01, 1, 1],
                   [-0.01, 1, 0],
                   [-0.01, 0, 0],
                   [-0.01, 0, 1]])*1.2
tris2=np.array([[0,1,2], [2,3,0]])
mesh2 = trimesh.Trimesh(points2, tris2)
for ii in range(7):
    mesh2 =mesh2.subdivide()

points3 = np.array([[1, 1, -0.01],
                   [1, 0, -0.01],
                   [0, 0, -0.01],
                   [0, 1, -0.01]])*1.2
mesh3 = trimesh.Trimesh(points3, tris2)
for ii in range(7):
    mesh3 =mesh3.subdivide()

points4 = np.array([[1, -0.01, 1],
                   [1, -0.01, 0],
                   [0, -0.01, 0],
                   [0,  -0.01, 1]])*1.2
mesh4 = trimesh.Trimesh(points4, tris2)
for ii in range(7):
    mesh4 =mesh4.subdivide()

# Difference vectors
RR2 = mesh2.vertices[:,None,None,:] - p_tris[None,:,:,:]
RR3 = mesh3.vertices[:,None,None,:] - p_tris[None,:,:,:]
RR4 = mesh4.vertices[:,None,None,:] - p_tris[None,:,:,:]

tn, ta = tri_normals_and_areas(points, tris)

#%% Plot potentials on the planes and the respective sources on the triangle
for ii, func in enumerate((triangle_potential_uniform,
                       triangle_potential_dipole_linear,
                       gamma0, omega)):

    mlab.figure(bgcolor=(1,1,1))
    # Plot shape and potential
    print(func)
    if ii==0:
        pot2 = func(RR2, tn)[:,0]
        pot3 = func(RR3, tn)[:,0]
        pot4 = func(RR4, tn)[:,0]
        mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(0.5,0.5,0.5),
                             opacity=0.5)
    if ii==1:oome
        pot2 = func(RR2, tn, ta)[:,0,1]
        pot3 = func(RR3, tn, ta)[:,0,1]
        pot4 = func(RR4, tn, ta)[:,0,1]
        mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(0.5,0.5,0.5),
                             opacity=0.2)
        meshq = mesh.copy()
        for jj in range(4):
            meshq =meshq.subdivide()
        u = np.zeros(meshq.vertices.shape)
        r = meshq.vertices
        r2 = np.zeros(r.shape +(3,))
        r2[:,1] = r
        r2[:,0] = mesh.vertices[0]
        r2[:,2] = mesh.vertices[2]
        u[:, 2] = np.linalg.det(r2)/np.linalg.det(mesh.vertices)
        mlab.quiver3d(*r.T, *u.T, colormap='viridis', mode='arrow')


    if ii==2:
        pot2 = func(RR2)[:,0,2]
        pot3 = func(RR3)[:,0,2]
        pot4 = func(RR4)[:,0,2]
        mlab.plot3d(*points[0:2].T, color=(0.5,0.5,0.5),
                             opacity=0.5, tube_radius=0.02)
    if ii==3:
        pot2 = func(RR2)[:,0]
        pot3 = func(RR3)[:,0]
        pot4 = func(RR4)[:,0]
        mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(0.5,0.5,0.5),
                             opacity=0.2)
        u = np.zeros(meshq.vertices.shape)
        u[:, 2] = 1
        r = meshq.vertices
        mlab.quiver3d(*r.T, *u.T, colormap='viridis', mode='arrow')


    M = max(max(abs(pot2)),max(abs(pot3)),max(abs(pot4)))
    for m, p in zip((mesh2, mesh3, mesh4), (pot2, pot3, pot4)):
        s = mlab.triangular_mesh(*m.vertices.T, m.faces, scalars=p,
                                 colormap='seismic', vmin = -M,
                                 vmax = M)
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.module_manager.scalar_lut_manager.number_of_colors = 32




