"""
Biplanar coil design
====================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes.

"""

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import trimesh


from bfieldtools.conductor import Conductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

import pkg_resources


# Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


# Load simple plane mesh that is centered on the origin
planemesh = trimesh.load(
    file_obj=pkg_resources.resource_filename(
        "bfieldtools", "example_meshes/10x10_plane_hires.obj"
    ),
    process=False,
)

planemesh.apply_scale(scaling_factor * 1.6)

# Specify coil plane geometry
center_offset = np.array([0, 0, 0]) * scaling_factor
standoff = np.array([0, 5, 0]) * scaling_factor

# Create coil plane pairs
coil_plus = trimesh.Trimesh(
    planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
)

coil_minus = trimesh.Trimesh(
    planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
)

joined_planes = coil_plus.union(coil_minus)

# Create mesh class object
coil = Conductor(mesh_obj=joined_planes, fix_normals=True, basis_name="suh", N_suh=100)

##############################################################
# Set up target and stray field points

# Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 0]) * scaling_factor

sidelength = 2 * scaling_factor
n = 8
xx = np.linspace(-sidelength / 2, sidelength / 2, n)
yy = np.linspace(-sidelength / 2, sidelength / 2, n)
zz = np.linspace(-sidelength / 2, sidelength / 2, n)
X, Y, Z = np.meshgrid(xx, yy, zz, indexing="ij")

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

target_points = np.array([x, y, z]).T

# Turn cube into sphere by rejecting points "in the corners"
target_points = (
    target_points[np.linalg.norm(target_points, axis=1) < sidelength / 2] + center
)


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
# Create bfield specifications used when optimizing the coil geometry

# The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 0] = target_field[:, 0] + 1

target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": 0.01,
    "target": target_field,
}
stray_spec = {
    "coupling": coil.B_coupling(stray_points),
    "abs_error": 0.001,
    "target": np.zeros((n_stray_points, 3)),
}

bfield_specification = [target_spec, stray_spec]

##############################################################
## Compute the optimal stream function, either using a numerical solver or regularized least squares
import mosek

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec, stray_spec],
    objective="minimum_inductive_energy",
    solver="MOSEK",
    solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
)


#############################################################
# Plot coil windings and target points

N_contours = 10

loops, loop_values = scalar_contour(coil.mesh, coil.s, N_contours=N_contours)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors="auto", figure=f)

B_target = coil.B_coupling(target_points) @ coil.s

mlab.quiver3d(*target_points.T, *B_target.T)


##############################################################
# Plot cross-section of magnetic field and magnetic potential of the discretized loops

from bfieldtools.line_magnetics import magnetic_field, scalar_potential

x = y = np.linspace(-12, 12, 250)
X, Y = np.meshgrid(x, y, indexing="ij")
points = np.zeros((X.flatten().shape[0], 3))
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

B = np.zeros_like(points)
U = np.zeros((points.shape[0],))
for loop_idx in range(len(loops)):
    B += magnetic_field(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)
    U += scalar_potential(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)

B = B.T[:2].reshape(2, x.shape[0], y.shape[0])
lw = np.sqrt(B[0] ** 2 + B[1] ** 2)
lw = 2 * lw / np.max(lw)

U = U.reshape(x.shape[0], y.shape[0])

plt.figure()

plt.pcolormesh(X, Y, U.T, cmap="seismic", shading="gouraud")
# plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
#           extent=(x.min(), x.max(), y.min(), y.max()))

seed_points = points[:, :2] * 0.3

plt.streamplot(
    x,
    y,
    B[1],
    B[0],
    density=2,
    linewidth=lw,
    color="k",
    integration_direction="both",
    start_points=seed_points,
)
plt.axis("equal")
plt.axis("off")
for loop in loops:
    plt.plot(loop[:, 1], loop[:, 0], "k", linewidth=4, alpha=0.1)

plt.tight_layout()


#######################################################
# Lets also do the same coil optimization using regularized least-squares


from bfieldtools.coil_optimize import optimize_lsq

coil.s2 = optimize_lsq(
    coil, [target_spec, stray_spec], objective="minimum_resistive_energy", reg=1e6
)


#############################################################
# Plot coil windings and target points

N_contours = 10

loops, loop_values = scalar_contour(coil.mesh, coil.s2, N_contours=N_contours)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors="auto", figure=f)

B_target = coil.B_coupling(target_points) @ coil.s2

mlab.quiver3d(*target_points.T, *B_target.T)


##############################################################
# Plot cross-section of magnetic field and magnetic potential of the discretized loops

from bfieldtools.line_magnetics import magnetic_field, scalar_potential

x = y = np.linspace(-12, 12, 250)
X, Y = np.meshgrid(x, y, indexing="ij")
points = np.zeros((X.flatten().shape[0], 3))
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

B = np.zeros_like(points)
U = np.zeros((points.shape[0],))
for loop_idx in range(len(loops)):
    B += magnetic_field(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)
    U += scalar_potential(np.vstack((loops[loop_idx], loops[loop_idx][0])), points)

B = B.T[:2].reshape(2, x.shape[0], y.shape[0])
lw = np.sqrt(B[0] ** 2 + B[1] ** 2)
lw = 2 * lw / np.max(lw)

U = U.reshape(x.shape[0], y.shape[0])

plt.figure()

plt.pcolormesh(X, Y, U.T, cmap="seismic", shading="gouraud")
# plt.imshow(U, vmin=-1.0, vmax=1.0, cmap='seismic', interpolation='bicubic',
#           extent=(x.min(), x.max(), y.min(), y.max()))

seed_points = points[:, :2] * 0.3

plt.streamplot(
    x,
    y,
    B[1],
    B[0],
    density=2,
    linewidth=lw,
    color="k",
    integration_direction="both",
    start_points=seed_points,
)
plt.axis("equal")
plt.axis("off")
for loop in loops:
    plt.plot(loop[:, 1], loop[:, 0], "k", linewidth=4, alpha=0.1)

plt.tight_layout()
