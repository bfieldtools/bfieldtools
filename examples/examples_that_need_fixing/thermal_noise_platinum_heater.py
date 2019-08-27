"""
Thin-film heater noise
=========================

This example computes the thermal magnetic noise produced by a platinum
thin-film heater geometry.
"""

import numpy as np
import matplotlib.pyplot as plt
import trimesh
from mayavi import mlab

from bfieldtools.thermal_noise import compute_current_modes, compute_dc_Bnoise

import pkg_resources

#%% Fix the simulation parameters and load the heater geometry

d = 1e-6
sigma = 1/(16.592 * 1e-6 * 1e-2) #Platinum @ 450 K
T = 450
kB = 1.38064852e-23
mu0 = 4*np.pi*1e-7


scale_factor = 1e-3



mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/thinspiral.obj'))

#Subdivide mesh for higher accuracy if needed
#mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)


# Center the mesh at the origin, apply scaling
mesh.apply_scale(scale_factor)

mesh.vertices[:, 2] = 0
mesh.vertices[:, 1] -= np.mean(mesh.vertices[:, 1])
mesh.vertices[:, 0] -= np.mean(mesh.vertices[:, 0])



#%% Compute the normalized thermal current modes, and thereafter compute the magnetic field noise
#   caused by the currents. Finally, visualize the result.   

vl = compute_current_modes(mesh)

Np = 30

zl = np.linspace(0.1, 5, Np) * scale_factor
fp = np.array((np.zeros(zl.shape), np.zeros(zl.shape)-0.001, zl)).T

B = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)

fig = plt.figure(figsize=(4, 3))

plt.semilogy(zl*1e3, np.linalg.norm(B, axis=1)*1e15, 'k')
plt.xlabel('Distance (mm)')
plt.ylabel('DC noise amplitude (fT/rHz)')

fig.tight_layout()

plt.grid()
fig.axes[0].spines['top'].set_visible(False)
fig.axes[0].spines['right'].set_visible(False)



mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
mlab.points3d(*fp.T)



Ngrid = 40
Ngrid_z = 5
xx = np.linspace(-8, 8, Ngrid) * scale_factor
yy = np.linspace(-8, 8, Ngrid) * scale_factor
zz = np.linspace(0, 3, Ngrid_z) * scale_factor
X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

grid_points = np.vstack((x, y, z)).T


mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
mlab.points3d(*grid_points.T)

B_grid = compute_dc_Bnoise(mesh, vl, grid_points, sigma, d, T)

B_grid_matrix = B_grid.reshape((Ngrid, Ngrid, Ngrid_z, 3))

B_grid_matrix_norm = np.linalg.norm(B_grid_matrix, axis=-1)

field = mlab.pipeline.vector_field(X, Y, Z, B_grid_matrix[:,:,:,0], B_grid_matrix[:,:,:,1], B_grid_matrix[:,:,:,2],
                              scalars=B_grid_matrix_norm, name='B-field')
#
#vectors = mlab.pipeline.vectors(field,
#                      scale_factor=(X[1, 0, 0] - X[0, 0, 0]),
#                      )
#
#
#vectors.glyph.mask_input_points = True
#vectors.glyph.mask_points.on_ratio = 5
#
#vcp = mlab.pipeline.vector_cut_plane(field)
#vcp.glyph.glyph.scale_factor=10*(X[1, 0, 0] - X[0, 0, 0])
## For prettier picture:
#vcp.implicit_plane.widget.enabled = True


iso = mlab.pipeline.iso_surface(field,
                                opacity=0.2,
                                colormap='viridis')

# A trick to make transparency look better: cull the front face
iso.actor.property.frontface_culling = False