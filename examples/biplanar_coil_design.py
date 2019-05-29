#%% Compact example of design of a biplanar coil


import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh


from bfieldtools.mesh_class import MeshWrapper
from bfieldtools.magnetic_field_mesh import compute_C
from bfieldtools.coil_optimize import optimize_streamfunctions

import pkg_resources


#Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1e2


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
coil = MeshWrapper(verts=joined_planes.vertices, tris=joined_planes.faces)

#%% Set up target and stray field points

#Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 0]) * scaling_factor

sidelength = 2 * scaling_factor
n = 6
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
stray_points_mesh = trimesh.creation.icosphere(subdivisions=2, radius=stray_radius)
stray_points = stray_points_mesh.vertices + center

n_stray_points = len(stray_points)



#%% Compute C matrices that are used to compute the generated magnetic field

coil.C = compute_C(coil.mesh, target_points, parallel=False)
coil.strayC = compute_C(coil.mesh, stray_points, parallel=False)

#%% Specify target field and run solver

#The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function
target_field = np.ones(target_points.shape[0], )


# The tolerance parameter will determine the spatial detail of the coil.
# Smaller tolerance means better but more intricate patterns. Too small values
# will not be solveable.
tolerance = 0.3

I, sol = optimize_streamfunctions(coil, target_field,
                         target_axis=0,
                         target_error={'on_axis':0.01, 'off_axis':0.01, 'stray':0.01},
                         laplacian_smooth=0,
                         tolerance=tolerance)


#%% Plot coil windings and target points

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
           size=(480, 480))
mlab.clf()

surface = mlab.pipeline.triangular_mesh_source(*coil.verts.T, coil.tris,scalars=I)

windings = mlab.pipeline.contour_surface(surface, contours=8)


B_target = np.vstack((coil.C[:, :, 0].dot(I),
                  coil.C[:, :, 1].dot(I),
                  coil.C[:, :, 2].dot(I))).T


mlab.quiver3d(*target_points.T, *B_target.T)


#%% Plot field falloff on two axes

plt.figure()

z1 = np.linspace(0, 30, 31) * scaling_factor

x1 = y1 = np.zeros_like(z1)

line1_points = np.vstack((x1, y1, z1)).T

line1_C = compute_C(coil.mesh, r=line1_points)

B_line1 = np.vstack((line1_C[:, :, 0].dot(I), line1_C[:, :, 1].dot(I), line1_C[:, :, 2].dot(I))).T

plt.semilogy(z1 / scaling_factor, np.linalg.norm(B_line1, axis=1)/np.mean(np.abs(target_field)), label='Z')

y2 = np.linspace(0, 30, 31) * scaling_factor

z2 = x2 = np.zeros_like(y2)

line2_points = np.vstack((x2, y2, z2)).T

line2_C = compute_C(coil.mesh, r=line2_points)

B_line2 = np.vstack((line2_C[:, :, 0].dot(I), line2_C[:, :, 1].dot(I), line2_C[:, :, 2].dot(I))).T

plt.semilogy(y2 / scaling_factor, np.linalg.norm(B_line2, axis=1)/np.mean(np.abs(target_field)), label='Y')
plt.ylabel('Field amplitude (target field units)')
plt.xlabel('Distance from origin')
plt.grid(True, which='minor', axis='y')
plt.grid(True, which='major', axis='y', color='k')
plt.grid(True, which='major', axis='x')

plt.legend()

plt.show()