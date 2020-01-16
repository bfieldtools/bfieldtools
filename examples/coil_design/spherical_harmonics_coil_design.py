'''
Spherical harmonics-generating coil design
==========================================

Example showing a basic biplanar coil producing a field profile defined by
spherical harmonics.

'''

#import sys
#path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
#if path in sys.path:
#    sys.path.insert(0, path)


import numpy as np
from mayavi import mlab
import trimesh

from bfieldtools.mesh_class import Conductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops


from bfieldtools.sphtools import compute_sphcoeffs_mesh, sphbasis, sphfittools
from bfieldtools import sphtools


import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


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

joined_planes = coil_plus.union(coil_minus)

#Create mesh class object
coil = Conductor(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)

lmax = 4
coil.C_alms, coil.C_blms = compute_sphcoeffs_mesh(coil.mesh, lmax=lmax)


#Radius of sphere of interest
Rmax = 1.0

lind = 0
coil.C_alms_norm = np.zeros_like(coil.C_alms)
for l in range(1,lmax+1):
    for m in range(-1*l,l+1):
        temp = (2*l**2 + l)*Rmax**(2*l-1)/(2*l-1)
        #coeffs2[lind] = coeffs[lind]**2*temp
        coil.C_alms_norm[lind] = coil.C_alms[lind]/temp**0.5
        lind += 1

target_alms = np.zeros((lmax * (lmax+2),))
target_blms = np.zeros((lmax * (lmax+2),))

target_blms[0] += 1


center = np.array([0, 0, 0]) * scaling_factor

sidelength = 2 * scaling_factor
n = 8
xx = np.linspace(-sidelength/2, sidelength/2, n)
yy = np.linspace(-sidelength/2, sidelength/2, n)
zz = np.linspace(-sidelength/2, sidelength/2, n)
X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

target_points = np.array([x, y, z]).T

#Turn cube into sphere by rejecting points "in the corners"
target_points = target_points[np.linalg.norm(target_points, axis=1) < sidelength/2]  + center




sph = sphbasis(4)
sphfield = sphtools.field(target_points, target_alms, target_blms, lmax)

target_field = sphfield/np.max(sphfield[:, 0])

target_field[:, 2] = 0

coil.plot_mesh()
mlab.quiver3d(*target_points.T, *sphfield.T)



##############################################################
# Create bfield specifications used when optimizing the coil geometry


target_spec = {'coupling':coil.C_blms, 'rel_error':0, 'abs_error':0.01, 'target':target_blms}


##############################################################
# Run QP solver
import mosek

coil.j, prob = optimize_streamfunctions(coil,
                                   [target_spec],
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )

B_target = coil.B_coupling(target_points) @ coil.j


lmax = 4
coil.C_alms, coil.C_blms = compute_sphcoeffs_mesh(coil.mesh, lmax=lmax)

Alms, Blms = coil.C_alms @ coil.j, coil.C_blms @ coil.j

Alms = np.zeros_like(Blms)
sphfield_target = sphtools.field(target_points, Alms, Blms, lmax)


coeffs, coeffs2, nrmse = sphfittools.fitSpectra(sph, np.repeat(target_points[:, :, None], 3, -1), B_target, lmax)



#############################################################
# Plot coil windings and target points

N_contours = 10

loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=N_contours)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f)

B_target = coil.B_coupling(target_points) @ coil.j

mlab.quiver3d(*target_points.T, *B_target.T)