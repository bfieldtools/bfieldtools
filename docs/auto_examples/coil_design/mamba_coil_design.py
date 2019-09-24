'''
MAMBA coil
==========

Compact example of a biplanar coil producing homogeneous field in a number of target
regions arranged in a grid.

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
planemesh = trimesh.load(file_obj=pkg_resources.resource_filename('bfieldtools', 'example_meshes/10x10_plane_hires.obj'), process=False)

planemesh.apply_scale(scaling_factor)

#planemesh.vertices, planemesh.faces = trimesh.remesh.subdivide(planemesh.vertices, planemesh.faces)


#Specify coil plane geometry
center_offset = np.array([0, 0, 0]) * scaling_factor
standoff = np.array([0, 1.5, 0]) * scaling_factor

#Create coil plane pairs
coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
                         planemesh.faces, process=False)

coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
                     planemesh.faces, process=False)

joined_planes = coil_plus.union(coil_minus)

#Create mesh class object
coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces, fix_normals=True)

###############################################################
# Set up target and stray field points. Here, the target points are on a planar
# 4x4 grid slightly smaller than the coil dimensions.

center = np.array([0, 0, 0]) * scaling_factor

sidelength = 0.5 * scaling_factor
n = 4

height = 0.1
n_height = 2
xx = np.linspace(-sidelength/2, sidelength/2, n)
yy = np.linspace(-height/2, height/2, n_height)
zz = np.linspace(-sidelength/2, sidelength/2, n)
X, Y, Z = np.meshgrid(xx, yy, zz, indexing='ij')

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

target_points = np.array([x, y, z]).T


grid_target_points = list()
target_field = list()

hori_offsets = [-3, -1, 1, 3]
vert_offsets = [-3, -1, 1, 3]

for i, offset_x in enumerate(hori_offsets):
    for j, offset_y in enumerate(vert_offsets):
        grid_target_points.append(target_points + np.array([offset_x, 0, offset_y]))
        target_field.append((i + j - 3) * np.ones((len(target_points),)))

target_points = np.asarray(grid_target_points).reshape((-1,3))
target_field = np.asarray(target_field).reshape((-1,))

target_field = np.array([np.zeros((len(target_field),)), target_field, np.zeros((len(target_field),))]).T

###############################################################
# Plot target points and mesh
scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))

mlab.quiver3d(*target_points.T, *target_field.T)
coil.plot_mesh()


###############################################################
# Compute C matrices that are used to compute the generated magnetic field, create field specification

coil.C = compute_C(coil.mesh, target_points)

target_spec = {'C':coil.C, 'rel_error':0.01, 'abs_error':0, 'target_field':target_field}

###############################################################
# Run QP solver


# The tolerance parameter will determine the spatial detail of the coil.
# Smaller tolerance means better but more intricate patterns. Too small values
# will not be solveable.
tolerance = 0.45

coil.I, coil.sol = optimize_streamfunctions(coil,
                                            [target_spec],
                                            objective='minimum_inductive_energy',
                                            tolerance=tolerance)


###############################################################
# Plot coil windings and target points

loops, loop_values= scalar_contour(coil.mesh, coil.I, N_contours=10)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors='auto', figure=f, tube_radius=0.01)

B_target = coil.C.transpose([0, 2, 1]) @ coil.I

mlab.quiver3d(*target_points.T, *B_target.T)

f.scene.isometric_view()
