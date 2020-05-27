"""
MAMBA coil
==========

Compact example of a biplanar coil producing homogeneous field in a number of target
regions arranged in a grid. Meant to demonstrate the flexibility in target choice, inspired by the 
technique "multiple-acquisition micro B(0) array" (MAMBA) technique, see https://doi.org/10.1002/mrm.10464

"""


import numpy as np
from mayavi import mlab
import trimesh


from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

from bfieldtools.utils import combine_meshes, load_example_mesh


# Load simple plane mesh that is centered on the origin
planemesh = load_example_mesh("10x10_plane_hires")

# Specify coil plane geometry
center_offset = np.array([0, 0, 0])
standoff = np.array([0, 1.5, 0])

# Create coil plane pairs
coil_plus = trimesh.Trimesh(
    planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
)

coil_minus = trimesh.Trimesh(
    planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
)

joined_planes = combine_meshes((coil_plus, coil_minus))

# Create mesh class object
coil = MeshConductor(mesh_obj=joined_planes, fix_normals=True, basis_name="inner")

#%%
# Set up target and stray field points. Here, the target points are on a planar
# 4x4 grid slightly smaller than the coil dimensions.

center = np.array([0, 0, 0])

sidelength = 0.5
n = 4

height = 0.1
n_height = 2
xx = np.linspace(-sidelength / 2, sidelength / 2, n)
yy = np.linspace(-height / 2, height / 2, n_height)
zz = np.linspace(-sidelength / 2, sidelength / 2, n)
X, Y, Z = np.meshgrid(xx, yy, zz, indexing="ij")

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

target_points = np.asarray(grid_target_points).reshape((-1, 3))
target_field = np.asarray(target_field).reshape((-1,))

target_field = np.array(
    [np.zeros((len(target_field),)), target_field, np.zeros((len(target_field),))]
).T


target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 1] += 0.1
target_abs_error[:, 0::2] += 0.1

#%%
# Plot target points and mesh
coil.plot_mesh(opacity=0.1)
mlab.quiver3d(*target_points.T, *target_field.T)


#%%
# Compute coupling matrix that is used to compute the generated magnetic field, create field specification


target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": target_abs_error,
    "target": target_field,
}

#%%
# Run QP solver, plot result

import mosek

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_inductive_energy",
    solver="MOSEK",
    solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
)


coil.s.plot()

coil.s.discretize(N_contours=10).plot_loops()
