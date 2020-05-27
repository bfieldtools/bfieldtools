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


from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.viz import plot_cross_section
from bfieldtools.utils import combine_meshes, load_example_mesh


# Load simple plane mesh that is centered on the origin
planemesh = load_example_mesh("10x10_plane_hires")

# Specify coil plane geometry
center_offset = np.array([0, 0, 0])
standoff = np.array([0, 3, 0])

# Create coil plane pairs
coil_plus = trimesh.Trimesh(
    planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
)

coil_minus = trimesh.Trimesh(
    planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
)

joined_planes = combine_meshes((coil_plus, coil_minus))

# Create mesh class object
coil = MeshConductor(
    mesh_obj=joined_planes, fix_normals=True, basis_name="suh", N_suh=100
)

#%%
# Set up target and stray field points

# Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 0])

sidelength = 1.5
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


# Here, the stray field points are on a spherical surface
stray_radius = 20
stray_points_mesh = trimesh.creation.icosphere(subdivisions=3, radius=stray_radius)
stray_points = stray_points_mesh.vertices + center

n_stray_points = len(stray_points)


#%%
# Create bfield specifications used when optimizing the coil geometry

# The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 0] += 1

target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": 0.01,
    "target": target_field,
}
stray_spec = {
    "coupling": coil.B_coupling(stray_points),
    "abs_error": 0.01,
    "target": np.zeros((n_stray_points, 3)),
}

bfield_specification = [target_spec, stray_spec]

#%%
## Compute the optimal stream function, either using a numerical solver or regularized least squares

import mosek

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec, stray_spec],
    objective="minimum_ohmic_power",
    solver="MOSEK",
    solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
)


#%%
# Plot the optimized stream function, then discretize it and plot coil windings and the resultant magnetic field

coil.s.plot()

loops = coil.s.discretize(N_contours=10)

loops.plot_loops()

B_target = loops.magnetic_field(target_points)
mlab.quiver3d(*target_points.T, *B_target.T)


#%%
# Lets also do the same coil optimization using regularized least-squares.
# Now we can't specify inequality constraints (e.g. use error margins in the specification).


from bfieldtools.coil_optimize import optimize_lsq

coil.s2 = optimize_lsq(
    coil, [target_spec, stray_spec], objective="minimum_ohmic_power", reg=1e6
)


#%%
# Plot the optimized stream function, then discretize it and plot coil windings and the resultant magnetic field

coil.s2.plot()

loops2 = coil.s2.discretize(N_contours=10)

loops2.plot_loops()

B_target = loops2.magnetic_field(target_points)
mlab.quiver3d(*target_points.T, *B_target.T)


#%%
# Plot cross-section of magnetic field and magnetic potential of the discretized loops

x = y = np.linspace(-12, 12, 250)
X, Y = np.meshgrid(x, y, indexing="ij")


points = np.zeros((X.flatten().shape[0], 3))
points[:, 0] = X.flatten()
points[:, 1] = Y.flatten()

B = loops2.magnetic_field(points)
U = loops2.scalar_potential(points)

U = U.reshape(x.shape[0], y.shape[0])
B = B.T[:2].reshape(2, x.shape[0], y.shape[0])

lw = np.sqrt(B[0] ** 2 + B[1] ** 2)

lw = 2 * lw / np.max(lw)

plot_cross_section(X, Y, U, log=False, contours=False)

seed_points = points[:, :2] * 0.3

plt.streamplot(
    x,
    y,
    B[0],
    B[1],
    density=2,
    linewidth=lw,
    color="k",
    integration_direction="both",
    start_points=seed_points,
)
plt.axis("equal")
plt.axis("off")

plt.plot([-5, 5], [-3, -3], "k", linewidth=3, alpha=1)
plt.plot([-5, 5], [3, 3], "k", linewidth=3, alpha=1)

plt.tight_layout()
