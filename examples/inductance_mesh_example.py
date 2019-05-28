"""
Example for calculating independent current modes in a wall with a door
"""

import numpy as np
from mayavi import mlab
from matplotlib.tri import Triangulation
from scipy.linalg import eigh


from bfieldtools.laplacian_mesh import laplacian_matrix
from bfieldtools.mutual_inductance_mesh import self_inductance_matrix
from bfieldtools.utils import tri_normals_and_areas, dual_areas

xx = np.linspace(0, 1, 50)
X,Y = np.meshgrid(xx,xx,indexing='ij')
x = X.ravel()
y = Y.ravel()
z = np.zeros_like(x)
print('Triangulating mesh')
tt = Triangulation(x,y)

verts = np.array([x,y,z]).T
tris = tt.triangles

# Determine inner triangles, this is a mess currently
# Could be done in a systematic manner
centers = np.mean([x[tris], y[tris]], axis=-1)
trimask = ~((centers[1] < 0.7)*(abs(centers[0]-0.2) < 0.1))
tris2 = tris[trimask]
tt2 = Triangulation(x,y, tris2)
boundary_tris = np.array([-1 in n for n in tt2.neighbors])
all_inds = np.unique(tris.flatten())
boundary_inds0 = np.unique(tris2[boundary_tris]).flatten()
boundary_inds = [b for b in boundary_inds0 if np.sum(tris2==b)<=4]
boundary_inds = np.array(list(boundary_inds) +
                  list(np.nonzero((abs(y-0.7)<0.01)*(abs(x-0.3)<0.01))[0]))
#    plt.triplot(x, y, tris)
#    plt.plot(x[boundary_inds], y[boundary_inds], 'r*')

print('Calculating triangle stuff')
n, a = tri_normals_and_areas(verts, tris)
da = dual_areas(tris, a)

boundary_all = (np.isclose(x,0))+(np.isclose(y,0))+(np.isclose(x,1))+(np.isclose(y,1)) > 0

inner_inds = np.setdiff1d(all_inds,boundary_inds)
inner_inds = np.setdiff1d(inner_inds,np.nonzero(boundary_all)[0])


M = self_inductance_matrix(verts, tris)
M=0.5*(M+M.T)
Min = M[inner_inds[None,:], inner_inds[:,None]]
print('Calculating modes')
L = laplacian_matrix(verts, tris)
L = np.array(L.todense())
w,v = eigh(-L[inner_inds[None,:], inner_inds[:,None]], Min)

#%% Plot eigenmodes of surface currents on thin wall
mlab.figure()
scalars = np.zeros(x.shape)
Nmodes = 16
limit = np.max(abs(v[:,0]))
for ii in range(Nmodes):
    n = int(np.sqrt(Nmodes))
    i = ii % n
    j = int(ii/n)
    print(i,j)
    x = verts[:,0] + i*1.1
    y = verts[:,1] + j*1.1
    z = verts[:,2]
    scalars[inner_inds] = v[:,ii]
    #        scalars[inner] = v[:,4] +v[:,5]
    s=mlab.triangular_mesh(x,y,z, tris, scalars=scalars) #M[:,70])
    s.module_manager.scalar_lut_manager.number_of_colors = 16
    s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
    s.actor.mapper.interpolate_scalars_before_mapping = True
