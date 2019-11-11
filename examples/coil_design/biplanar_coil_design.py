'''
Biplanar coil design
====================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes.

'''

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh

from bfieldtools.mesh_class import MeshWrapper, CouplingMatrix
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


#Load simple plane mesh that is centered on the origin
planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane.obj'), process=False)

planemesh.apply_scale(scaling_factor*1.6)

#Specify coil plane geometry
center_offset = np.array([0, 0, 0]) * scaling_factor
standoff = np.array([0, 5, 0]) * scaling_factor

#Create coil plane pairs
coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                         planemesh.faces, process=False)

coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                     planemesh.faces, process=False)

joined_planes = coil_plus.union(coil_minus)

#Create mesh class object
coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)

##############################################################
# Set up target and stray field points

#Here, the target points are on a volumetric grid within a sphere

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



#    #Here, the stray field points are on a spherical surface
stray_radius = 20 * scaling_factor
#    stray_length = 20 * scaling_factor
#
#    stray_points = cylinder_points(radius=stray_radius,
#                                   length = stray_length,
#                                   nlength = 5,
#                                   nalpha = 30,
#                                   orientation=np.array([1, 0, 0]))
#
stray_points_mesh = trimesh.creation.icosphere(subdivisions=3, radius=stray_radius)
stray_points = stray_points_mesh.vertices + center

n_stray_points = len(stray_points)



##############################################################
# Compute C matrices that are used to compute the generated magnetic field

coil.C = CouplingMatrix(coil, compute_C)

##############################################################
# Create bfield specifications used when optimizing the coil geometry

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 0] = target_field[:, 0] + 1

target_rel_error = np.zeros_like(target_field)
target_rel_error[:, 0] += 0.01

target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 0] += 0.001
target_abs_error[:, 1:3] += 0.005

target_spec = {'C':coil.C(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target_field':target_field}
stray_spec = {'C':coil.C(stray_points), 'abs_error':0.01, 'rel_error':0, 'target_field':np.zeros((n_stray_points, 3))}

bfield_specification = [target_spec, stray_spec]

##############################################################
# Run QP solver
import mosek

coil.I, prob = optimize_streamfunctions(coil,
                                   [target_spec, stray_spec],
                                   objective='minimum_inductive_energy',
                                   solver='CVXOPT',
                                   solver_opts={}
                                   #{'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )

#############################################################
# Plot coil windings and target points

N_contours = 10

loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=N_contours)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f)

B_target = coil.C(target_points) @ coil.I

mlab.quiver3d(*target_points.T, *B_target.T)