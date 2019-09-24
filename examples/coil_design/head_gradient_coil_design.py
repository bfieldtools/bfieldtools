'''
Head gradient coil
==================

Example showing a gradient coil designed on the surface of a MEG system helmet
'''


import numpy as np
from mayavi import mlab
import trimesh


from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


#Load simple plane mesh that is centered on the origin
helmetmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools',
                                                                   'example_meshes/meg_helmet.obj'),
                          process=False)

#planemesh.apply_scale(scaling_factor)
#
##Specify coil plane geometry
#center_offset = np.array([0, 0, 0]) * scaling_factor
#standoff = np.array([0, 4, 0]) * scaling_factor
#
##Create coil plane pairs
#coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
#                         planemesh.faces, process=False)
##
#coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
#                     planemesh.faces, process=False)

#joined_planes = coil_plus.union(coil_minus)

#Create mesh class object
coil = MeshWrapper(verts=helmetmesh.vertices, tris=helmetmesh.faces, fix_normals=True)

###############################################################
#Set up target and stray field points.
#Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 0.04]) * scaling_factor

sidelength = 0.1 * scaling_factor
n = 12
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


###############################################################
#Compute C matrices that are used to compute the generated magnetic field

coil.C = compute_C(coil.mesh, target_points)


###############################################################
#Specify target field and run solver

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function
target_field = np.ones(target_points.shape[0], ) * target_points[:,0]



target_field = np.zeros(target_points.shape)
target_field[:, 0] = target_field[:, 0] + 1 * target_points[:,0]/np.max(target_points[:,0])

target_spec = {'C':coil.C, 'rel_error':0.1, 'abs_error':0, 'target_field':target_field}

# The tolerance parameter will determine the spatial detail of the coil.
# Smaller tolerance means better but more intricate patterns. Too small values
# will not be solveable.
tolerance = 5

coil.I, coil.sol = optimize_streamfunctions(coil,
                                            [target_spec],
                                            objective='minimum_inductive_energy',
                                            tolerance=tolerance)


###############################################################
#Plot coil windings and magnetic field in target points


loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=10)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.05/100)

B_target = coil.C.transpose([0, 2, 1]) @ coil.I

mlab.quiver3d(*target_points.T, *B_target.T)

f.scene.isometric_view()
