#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 10:15:23 2019

@author: makinea1
"""

import sys

path = "/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools"
if path in sys.path:
    sys.path.insert(0, path)


import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh
import pkg_resources

from bfieldtools.mesh_magnetics import magnetic_field_coupling
from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic
from bfieldtools.mesh_magnetics import scalar_potential_coupling
from bfieldtools.sphtools import compute_sphcoeffs_mesh, sphbasis
from bfieldtools.suhtools import suhbasis
from bfieldtools.integrals import triangle_potential_uniform


def single_layer_pot(mesh, points):
    """ Calculate single-layer potential for
        constant basis functions on triangles
    """
    R = points[:, None, None, :] - mesh.vertices[None, mesh.faces, :]

    u = triangle_potential_uniform(R, mesh.face_normals)

    return u


mesh = trimesh.load(
    file_obj=pkg_resources.resource_filename(
        "bfieldtools", "example_meshes/cube_fillet.stl"
    ),
    process=True,
)

mesh.vertices -= mesh.vertices.mean(axis=0)

# Points insided and outside the surface
mesh_in = mesh.copy()
mesh_in.vertices -= 0.005 * mesh_in.vertex_normals
mesh_out = mesh.copy()
mesh_out.vertices += 0.005 * mesh_in.vertex_normals

# Create some source pattern
bsuh = suhbasis(mesh, 20)
coeffs = np.zeros(20)
coeffs[2] = 1
coeffs[6] = 1
dd = bsuh.basis @ coeffs
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=dd)
# Csuh = magnetic_field_coupling_analytic(mesh, mesh_field.vertices) @ bsuh.basis

N = 50
x = np.linspace(-7, 7, N)
R = np.array(np.meshgrid(x, x, 0, indexing="ij")).reshape(3, N * N)

Ud_eval = scalar_potential_coupling(mesh, R.T)
ud = (Ud_eval @ dd).reshape(N, N)

#%% Solve for potential that is zero inside
points = mesh_in.vertices[mesh_in.faces].mean(axis=1)
uin_dd = scalar_potential_coupling(mesh, points)
uin_ss = single_layer_pot(mesh, points)
ss = -np.linalg.solve(uin_ss, uin_dd @ dd)

Us_eval = single_layer_pot(mesh, R.T)
us1 = (Us_eval @ ss).reshape(N, N)

#%% Solve for potential that is zero outside
points = mesh_out.vertices[mesh_out.faces].mean(axis=1)
uout_dd = scalar_potential_coupling(mesh, points)
uout_ss = single_layer_pot(mesh, points)
ss_out = -np.linalg.solve(uout_ss, uout_dd @ dd)

us2 = (Us_eval @ ss_out).reshape(N, N)

#%%
u1 = ud + us1  # Double and single layer that cancel each other inside
u2 = ud + us2  # Double and single layer that cancel each other outside

# Now utot2 - utot1 = us2 - us1 is a single layer potential
# That is equivalent to -utot1 outside, and utot2 inside
fig, axes = plt.subplots(1, 3)
plt.sca(axes[0])
plt.imshow(-u1)
plt.sca(axes[1])
plt.imshow(u2)
plt.sca(axes[2])
plt.imshow(us2 - us1)


#%% Solve for double-layer density dd that cancels a single-layer potential
points = mesh_in.vertices
uin_dd = scalar_potential_coupling(mesh, points)
uin_ss = single_layer_pot(mesh, points)
# Solve for dd1
dd1 = -np.linalg.solve(uin_dd, uin_ss @ ss_out)
ud1 = (Ud_eval @ dd1).reshape(N, N)

#%%
v1 = ud1 + us2  # Double and single layer that cancel each other inside
v2 = ud + us2  # Double and single layer that cancel each other outside

# Now utot2 - utot1 = us2 - us1 is a single layer potential
# That is equivalent to utot1 outside, and utot2 inside
fig, axes = plt.subplots(1, 3)
plt.sca(axes[0])
plt.imshow(-v1)
plt.sca(axes[1])
plt.imshow(v2)
plt.sca(axes[2])
plt.imshow(ud - ud1)
