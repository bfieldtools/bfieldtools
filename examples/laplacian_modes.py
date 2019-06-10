"""
Example for calculating and plotting eigenmodes of Laplacian matrix of surface mesh
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from mayavi import mlab
import trimesh

from bfieldtools import utils
from bfieldtools.laplacian_mesh import laplacian_matrix, mass_matrix

import pkg_resources

mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane.obj'))

mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

L = laplacian_matrix(mesh)
M = mass_matrix(mesh)

u, v = eigh(-L.todense()[inner_verts][:,inner_verts], M.todense()[inner_verts][:,inner_verts])

plt.plot(u)
scalars = np.zeros(L.shape[0])
scalars[inner_verts] = v[:, 0]
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=scalars)

Nmodes = 16**2
limit = np.max(abs(v[:,0]))

verts= mesh.vertices
tris = mesh.faces

for ii in range(Nmodes):
	n = int(np.sqrt(Nmodes))
	i = ii % n
	j = int(ii/n)
	print(i,j)
	x = verts[:,0] + i*12
	y = verts[:,1]
	z = verts[:,2] + j*12
	scalars[inner_verts] = v[:,ii]
	#        scalars[inner] = v[:,4] +v[:,5]
	s=mlab.triangular_mesh(x,y,z, tris, scalars=scalars) #M[:,70])
	s.module_manager.scalar_lut_manager.number_of_colors = 256
	s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
	s.actor.mapper.interpolate_scalars_before_mapping = True



# If you want, animate shift between modes
#    s = mlab.triangular_mesh(*verts.T, tris, scalars=scalars)
#    s.actor.mapper.interpolate_scalars_before_mapping = True
#
#    @mlab.animate
#    def anim():
#        for i in range(100000):
#            scalars = np.zeros(L.shape[0])
#
#            prev = int(i / 10)
#            post = prev + 1
#
#            trans = (i % 10) / 10
#
#            scalars[inner_verts] = (1 - trans ) * v[:, prev]  + trans * v[:, post]
#            s.mlab_source.scalars = np.asarray(scalars, 'd')
#            yield
#
#    anim()
#    mlab.show()