'''
Shielded cylindrical coil
=========================
Compact example of design of a shielded cylindrical coil, for which the transient field
due to inductive interaction with the shield is minimized
'''


import numpy as np
from mayavi import mlab
import trimesh


from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.mutual_inductance_mesh import mutual_inductance_matrix

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

coil.C = compute_C(coil.mesh, target_points)
shield.C = compute_C(shield.mesh, target_points)

mutual_inductance = mutual_inductance_matrix(coil.mesh, shield.mesh)

# Take into account the field produced by currents induced into the shield
# NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

shield.coupling = -np.linalg.pinv(shield.inductance) @ mutual_inductance.T
secondary_C = (shield.C.transpose((0,2,1)) @ shield.coupling).transpose((0,2,1))

###############################################################
# Create bfield specifications used when optimizing the coil geometry

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function
target_field = np.zeros(target_points.shape)
target_field[:, 1] = target_field[:, 1] + 1 # Homogeneous Z-field

target_spec = {'C':coil.C, 'rel_error':0.01, 'abs_error':0, 'target_field':target_field}


induction_spec = {'C':secondary_C, 'abs_error':0.1, 'rel_error':0, 'target_field':np.zeros(target_field.shape)}

###############################################################
# Run QP solver

# The tolerance parameter will determine the spatial detail of the coil.
# Smaller tolerance means better but more intricate patterns. Too small values
# will not be solveable.
tolerance = 0.5

coil.I, coil.sol = optimize_streamfunctions(coil,
                                            [target_spec, induction_spec],
                                            objective='minimum_inductive_energy',
                                            tolerance=tolerance)

shield.induced_I = shield.coupling @ coil.I


###############################################################
# Plot coil windings and target points

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

surface = mlab.pipeline.triangular_mesh_source(*coil.mesh.vertices.T, coil.mesh.faces,scalars=coil.I)

windings = mlab.pipeline.contour_surface(surface, contours=20)
windings.module_manager.scalar_lut_manager.number_of_colors = 2 #Color windings according to current direction
windings.module_manager.scalar_lut_manager.reverse_lut = True #Flip LUT for the colors to correspond to RdBu colormap

shield_surface = mlab.pipeline.triangular_mesh_source(*shield.mesh.vertices.T, shield.mesh.faces,scalars=shield.induced_I)

shield_surface_render = mlab.pipeline.surface(shield_surface, colormap='RdBu')

shield_surface_render.actor.property.frontface_culling = True

B_target = coil.C.transpose([0, 2, 1]) @ coil.I

mlab.quiver3d(*target_points.T, *B_target.T)

mlab.title('Coils which minimize the transient effects of conductive shield')


###############################################################
# For comparison, let's see how the coils look when we ignore the conducting shield


# The tolerance parameter will determine the spatial detail of the coil.
# Smaller tolerance means better but more intricate patterns. Too small values
# will not be solveable.
tolerance = 0.5

coil.unshielded_I, coil.unshielded_sol = optimize_streamfunctions(coil,
                                            [target_spec],
                                            laplacian_smooth=0,
                                            tolerance=tolerance)

shield.unshielded_induced_I = shield.coupling @ coil.unshielded_I

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
                size=(800, 800))
mlab.clf()

surface = mlab.pipeline.triangular_mesh_source(*coil.mesh.vertices.T, coil.mesh.faces,scalars=coil.unshielded_I)

windings = mlab.pipeline.contour_surface(surface, contours=20)
windings.module_manager.scalar_lut_manager.number_of_colors = 2 #Color windings according to current direction
windings.module_manager.scalar_lut_manager.reverse_lut = True #Flip LUT for the colors to correspond to RdBu colormap

shield_surface = mlab.pipeline.triangular_mesh_source(*shield.mesh.vertices.T, shield.mesh.faces,scalars=shield.unshielded_induced_I)

shield_surface_render = mlab.pipeline.surface(shield_surface, colormap='RdBu')

shield_surface_render.actor.property.frontface_culling = True

B_target = coil.C.transpose([0, 2, 1]) @ coil.unshielded_I

mlab.quiver3d(*target_points.T, *B_target.T)

mlab.title('Coils which ignore the conductive shield')

