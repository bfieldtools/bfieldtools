"""
Short example showing that the numerical field computation gives the same result as analytical computation for simple unit disc.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from mayavi import mlab
import trimesh

import bfieldtools
from bfieldtools import utils
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.laplacian_mesh import laplacian_matrix, mass_matrix

#    mesh = trimesh.load('./example_meshes/10x10_plane_hires.obj')
mesh = trimesh.load(os.path.join(bfieldtools.__file__,  './example_meshes/unit_disc.stl'))

mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
#    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
#    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

Np = 20
dz = 0.2

z = np.linspace(0.1, 10, Np)
fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T
mlab.figure()
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
mlab.points3d(*fp.T)

C = compute_C(mesh, fp)

k = 1

#    plt.figure()

r = np.linalg.norm(mesh.vertices,axis=1)
sf = -1*k*r

mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=sf)

B = np.zeros(fp.shape)
B[:,0] = C[:,:,0].dot(sf)
B[:,1] = C[:,:,1].dot(sf)
B[:,2] = C[:,:,2].dot(sf)


R = 1
mu0 = 4*np.pi*1e-7
Ban = 0.5*mu0*k*(np.log(np.sqrt(R**2+z**2)+R) - R/np.sqrt(R**2+z**2) - np.log(z))


plt.figure()
plt.plot(z, B[:,2])
plt.plot(z, Ban)


plt.figure()
plt.plot(z, B[:,2]/np.max(B[:,2]))
plt.plot(z, Ban/np.max(Ban))

plt.figure()
plt.plot(z, np.abs((B[:,2]-Ban))/np.abs(Ban)*100)