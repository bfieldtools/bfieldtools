"""
Example for calculating independent current modes on planar structure
"""

import numpy as np
from mayavi import mlab
import trimesh
from scipy.linalg import eigh
import pkg_resources


from bfieldtools.laplacian_mesh import laplacian_matrix
from bfieldtools.mutual_inductance_mesh import self_inductance_matrix
from bfieldtools.utils import tri_normals_and_areas, dual_areas, find_mesh_boundaries

# Load simple plane mesh that is centered on the origin
file_obj = file_obj = pkg_resources.resource_filename(
    "bfieldtools", "example_meshes/10x10_plane_hires.obj"
)

planemesh = trimesh.load(file_obj, process=False)

verts = planemesh.vertices
tris = planemesh.faces

print("Calculating triangle stuff")
n, a = tri_normals_and_areas(verts, tris)
da = dual_areas(tris, a)


M = self_inductance_matrix(planemesh)

boundary_verts, inner_verts, boundary_tris, inner_tris = find_mesh_boundaries(
    planemesh.vertices, planemesh.faces, planemesh.edges
)


M = 0.5 * (M + M.T)
Min = M[inner_verts[None, :], inner_verts[:, None]]
print("Calculating modes")
L = laplacian_matrix(planemesh)
L = np.array(L.todense())
w, v = eigh(-L[inner_verts[None, :], inner_verts[:, None]], Min)

#%% Plot eigenmodes of surface currents on thin wall
mlab.figure()
scalars = np.zeros(M.shape[0])
Nmodes = 16
limit = np.max(abs(v[:, 0]))
for ii in range(Nmodes):
    n = int(np.sqrt(Nmodes))
    i = ii % n
    j = int(ii / n)
    print(i, j)
    x = verts[:, 0] + i * 15
    y = verts[:, 1]
    z = verts[:, 2] + j * 15
    scalars[inner_verts] = v[:, ii]
    #        scalars[inner] = v[:,4] +v[:,5]
    s = mlab.triangular_mesh(x, y, z, tris, scalars=scalars)  # M[:,70])
    s.module_manager.scalar_lut_manager.number_of_colors = 16
    s.module_manager.scalar_lut_manager.data_range = np.array([-limit, limit])
    s.actor.mapper.interpolate_scalars_before_mapping = True
