# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 09:37:19 2019

@author: Antti
"""

import numpy as np
import trimesh
from timeit import timeit
from mayavi import mlab
import matplotlib.pyplot as plt

import pkg_resources
import sys
#path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
#if path not in sys.path:
#    sys.path.insert(0,path)

from bfieldtools.laplacian_mesh import laplacian_matrix, mass_matrix
from bfieldtools.mesh_class import MeshWrapper


#Load simple plane mesh that is centered on the origin
file_obj = pkg_resources.resource_filename('bfieldtools',
                    'example_meshes/plane_w_holes.stl')
coilmesh = trimesh.load(file_obj, process=True)
coil = MeshWrapper(mesh_obj = coilmesh)

mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces)
mlab.points3d([-6,6,0], [0,0,0], [0,0,0])
mlab.show()

#%%
i1 = np.sum((coilmesh.vertices[coil.boundary_verts] - np.array([6,0,0]))**2, axis=1) < 3**2
i2 = np.sum((coilmesh.vertices[coil.boundary_verts] - np.array([0,0,0]))**2, axis=1) < 3**2
i3 = np.sum((coilmesh.vertices[coil.boundary_verts] - np.array([-6,0,0]))**2, axis=1) < 3**2

b1 = coil.boundary_verts[i1]
b2 = coil.boundary_verts[i2]
b3 = coil.boundary_verts[i3]
b4 = coil.boundary_verts[np.logical_not(i1+i2+i3)]

#%%
L = laplacian_matrix(coilmesh)
M = mass_matrix(coilmesh)
Linner = L[coil.inner_verts, :][:, coil.inner_verts]
Minner = M[coil.inner_verts, :][:, coil.inner_verts]

L1 = np.array(np.sum(L[b1,:][:, coil.inner_verts], axis=0))
L2 = np.array(np.sum(L[b2,:][:, coil.inner_verts], axis=0))
L3 = np.array(np.sum(L[b3,:][:, coil.inner_verts], axis=0))

u1 = np.linalg.solve(Linner.toarray() , -L1[0])

scalars = np.zeros(L.shape[0])
scalars[coil.inner_verts] =u1
mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=scalars)
mlab.show()


#%%
L_holes = np.concatenate((Linner.toarray(), L1.T, L2.T, L3.T), axis=1)
L11 = np.concatenate((L1, np.array([[-L1.sum(), 0, 0]])), axis=1)
L21 = np.concatenate((L2, np.array([[0, -L2.sum(), 0]])), axis=1)
L31 = np.concatenate((L3, np.array([[0, 0, -L3.sum()]])), axis=1)
L_holes = np.concatenate((L_holes, L11, L21, L31) ,axis=0)

#%%
m = M.diagonal()
M_holes = np.diag(np.concatenate((Minner.diagonal(), np.array([np.sum(m[b1]),
                                 np.sum(m[b2]), np.sum(m[b3])]))))

#%%
from scipy.linalg import eigh
ee, vv = eigh(-L_holes, M_holes, eigvals=(0,10))

#%%
scalars = np.zeros(L.shape[0])
scalars[coil.inner_verts] = vv[:-3, 2]
mlab.triangular_mesh(*coilmesh.vertices.T, coilmesh.faces, scalars=scalars)
mlab.show()