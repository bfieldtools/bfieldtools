'''
Coil with minimal eddy currents
===============================
Compact example of design of a cylindrical coil surrounded by a RF shield, i.e. a conductive surface.
The effects of eddy currents due to inductive interaction with the shield is minimized
'''


import numpy as np
from mayavi import mlab
import trimesh


from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.mesh_properties import mutual_inductance_matrix
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


#Load example coil mesh that is centered on the origin
coilmesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/cylinder.stl'), process=True)

coilmesh.apply_scale(0.75)

coilmesh.vertices, coilmesh.faces = trimesh.remesh.subdivide(coilmesh.vertices, coilmesh.faces)

#Specify offset from origin
center_offset = np.array([0, 0, 0.75])

#Apply offset
coilmesh = trimesh.Trimesh(coilmesh.vertices + center_offset,
                            coilmesh.faces, process=False)

#Create mesh class object
coil = MeshWrapper(verts=coilmesh.vertices, tris=coilmesh.faces, fix_normals=True)

# Separate object for shield geometry
shield = MeshWrapper(mesh_file=pkg_resources.resource_filename('bfieldtools', 'example_meshes/cylinder.stl'), process=True, fix_normals=True)

###############################################################
# Set up target  points and plot geometry

#Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 3])

sidelength = 0.75 * scaling_factor
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


#Plot coil, shield and target points

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                size=(800, 800))

coil.plot_mesh()
shield.plot_mesh()
mlab.points3d(*target_points.T)




###############################################################
# Compute C matrices that are used to compute the generated magnetic field


mutual_inductance = mutual_inductance_matrix(coil.mesh, shield.mesh)

# Take into account the field produced by currents induced into the shield
# NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

shield.coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
secondary_C = shield.B_coupling(target_points) @ shield.coupling

###############################################################
# Create bfield specifications used when optimizing the coil geometry

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function
target_field = np.zeros(target_points.shape)
target_field[:, 1] = target_field[:, 1] + 1

target_rel_error = np.zeros_like(target_field)
target_rel_error[:, 1] += 0.01

target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 1] += 0.001
target_abs_error[:, 0::2] += 0.005

target_spec = {'coupling':coil.B_coupling(target_points), 'rel_error':target_rel_error, 'abs_error':target_abs_error, 'target':target_field}


induction_spec = {'coupling':secondary_C, 'abs_error':0.1, 'rel_error':0, 'target':np.zeros(target_field.shape)}

###############################################################
# Run QP solver

import mosek

coil.j, prob = optimize_streamfunctions(coil,
                                   [target_spec, induction_spec],
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )

shield.induced_j = shield.coupling @ coil.j


###############################################################
# Plot coil windings and target points


loops, loop_values= scalar_contour(coil.mesh, coil.j, N_contours=10)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.02)

B_target = coil.B_coupling(target_points) @ coil.j

mlab.quiver3d(*target_points.T, *B_target.T)


mlab.title('Coils which minimize the transient effects of conductive shield')


###############################################################
# For comparison, let's see how the coils look when we ignore the conducting shield


coil.unshielded_j, coil.unshielded_prob = optimize_streamfunctions(coil,
                                   [target_spec],
                                   objective='minimum_inductive_energy',
                                   solver='MOSEK',
                                   solver_opts={'mosek_params':{mosek.iparam.num_threads: 8}}
                                   )

shield.unshielded_induced_j = shield.coupling @ coil.unshielded_j

loops, loop_values= scalar_contour(coil.mesh, coil.unshielded_j, N_contours=10)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.02)

B_target_unshielded = coil.B_coupling(target_points) @ coil.unshielded_j

mlab.quiver3d(*target_points.T, *B_target_unshielded.T)

mlab.title('Coils which ignore the conductive shield')

