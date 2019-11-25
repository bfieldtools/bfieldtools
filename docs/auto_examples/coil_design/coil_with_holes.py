'''
Coil with interior holes
========================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes. The coil planes have holes in them,

'''

import numpy as np
from mayavi import mlab
import trimesh

from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


#Load simple plane mesh that is centered on the origin
planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/plane_w_holes.stl'), process=False)

angle=np.pi/2
rotation_matrix = np.array([[1, 0, 0, 0],
                            [0, np.cos(angle), -np.sin(angle), 0],
                            [0, np.sin(angle), np.cos(angle), 0],
                            [0, 0, 0, 1]
                              ])

planemesh.apply_transform(rotation_matrix)
planemesh.apply_scale(scaling_factor)

#Specify coil plane geometry
center_offset = np.array([0, 0, 0]) * scaling_factor
standoff = np.array([0, 20, 0]) * scaling_factor

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

sidelength = 10 * scaling_factor
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




####################################################################
# Let's find and separate the inner and outer boundaries of the coil mesh

inner_bounds = np.intersect1d(coil.boundary_verts, np.where(np.linalg.norm(coil.mesh.vertices[:,0::2], axis=1)< 0.015)[0])

centre_hole1 = np.intersect1d(np.intersect1d(coil.boundary_verts,
                                             np.where(np.linalg.norm(coil.mesh.vertices[:,0::2], axis=1)< 0.004)[0]),
                              np.where(coil.mesh.vertices[:,1] < 0)[0])

centre_hole2 = np.intersect1d(np.intersect1d(coil.boundary_verts,
                                             np.where(np.linalg.norm(coil.mesh.vertices[:,0::2], axis=1)< 0.004)[0]),
                              np.where(coil.mesh.vertices[:,1] > 0)[0])

left_hole1 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                           np.where(coil.mesh.vertices[:,0] < -0.004)[0]), inner_bounds),
                            np.where(coil.mesh.vertices[:,1] < 0)[0])

left_hole2 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                           np.where(coil.mesh.vertices[:,0] < -0.004)[0]), inner_bounds),
                            np.where(coil.mesh.vertices[:,1] > 0)[0])

right_hole1 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                           np.where(coil.mesh.vertices[:,0] > 0.004)[0]), inner_bounds),
                            np.where(coil.mesh.vertices[:,1] < 0)[0])

right_hole2 = np.intersect1d(np.intersect1d(np.intersect1d(coil.boundary_verts,
                                           np.where(coil.mesh.vertices[:,0] > 0.004)[0]), inner_bounds),
                            np.where(coil.mesh.vertices[:,1] > 0)[0])

outer_bounds = np.setdiff1d(coil.boundary_verts, inner_bounds)


graph = trimesh.graph.vertex_adjacency_graph(coil.mesh)

zero_eq_indices = outer_bounds
iso_eq_indices = [left_hole1, centre_hole1, right_hole1, left_hole2, centre_hole2, right_hole2]

boundary_constraints = {'zero_eq_indices':zero_eq_indices , 'iso_eq_indices': iso_eq_indices}

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

target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}

bfield_specification = [target_spec]

##############################################################
# Run QP solver
import mosek

coil.j, prob = optimize_streamfunctions(coil,
                                   bfield_specification,
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}},
                                   boundary_constraints=boundary_constraints
                                   )

#############################################################
# Plot coil windings and target points

N_contours = 20

loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=N_contours)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.1)

B_target = coil.B_coupling(target_points) @ coil.j

mlab.quiver3d(*target_points.T, *B_target.T)