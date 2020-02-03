#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:00:56 2020

@author: makinea1
"""

import sys
path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
if path not in sys.path:
    sys.path.insert(0, path)


import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh

from bfieldtools.mesh_class import Conductor, StreamFunction
from bfieldtools.mesh_magnetics import magnetic_field_coupling as compute_C
from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic as compute_C_analytic
from bfieldtools.mesh_magnetics import scalar_potential_coupling as compute_U
from bfieldtools.mesh_properties import mutual_inductance_matrix
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops
from bfieldtools.sphtools import compute_sphcoeffs_mesh, ylm, cartesian2spherical
from bfieldtools.utils import dual_areas

from trimesh.creation import icosphere

# Volume of interest
r=0.2
target = icosphere(4, r)


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

coil1 = Conductor(mesh_obj=mesh1, basis_name = 'inner', N_sph = 4)

mlab.triangular_mesh(*mesh1.vertices.T, mesh1.faces)
mlab.triangular_mesh(*target.vertices.T, target.faces)

#%%
B = coil1.B_coupling(target.vertices)
U = coil1.U_coupling(target.vertices)
A = dual_areas(target.faces, target.area_faces)

Bn = np.sum(B*target.vertex_normals[:,:,None], axis=1)

# Field energy matrix
E = (Bn*A[:,None]).T @ U

from scipy.linalg import eigh
#evals, evecs = eigh(E, coil1.mass, eigvals=(0, 20))
evals, evecs = eigh(E, coil1.inductance, eigvals=(0, 20))

#%%
s = StreamFunction(evecs[:,5], coil1)
mlab.figure()
s.plot(True)
mlab.triangular_mesh(*target.vertices.T, target.faces, scalars = Bn @ s)
vectors=mlab.quiver3d(*target.vertices.T, *(B @ s).T, mode='arrow', color=(0,0,0))
vectors.glyph.glyph_source.glyph_position = 'center'

#%%
sph_coords =cartesian2spherical(target.vertices)
theta = sph_coords[:,1]
phi = sph_coords[:,2]
p = []
for l in range(1,5):
    for m in range(-l, l+1):
        p.append(ylm(l, m, theta, phi))
p = np.array(p)
#%%
#mlab.triangular_mesh(*target.vertices.T, target.faces, scalars = p[0])

#%% "Mutual energy" with sph components
innerprods = p @ ((A[:,None]*Bn) @ evecs)
plt.matshow(innerprods/np.linalg.norm(innerprods, axis=0, keepdims=True))
plt.xlabel('Eigen component')
plt.ylabel('SPH')


