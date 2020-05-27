"""
Coil with interior holes
========================

Example showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes. The coil planes have holes in them,

"""

import numpy as np
import trimesh

from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.utils import combine_meshes, load_example_mesh


# Load simple plane mesh that is centered on the origin
planemesh = load_example_mesh("plane_w_holes")

angle = np.pi / 2
rotation_matrix = np.array(
    [
        [1, 0, 0, 0],
        [0, np.cos(angle), -np.sin(angle), 0],
        [0, np.sin(angle), np.cos(angle), 0],
        [0, 0, 0, 1],
    ]
)

planemesh.apply_transform(rotation_matrix)

# Specify coil plane geometry
center_offset = np.array([0, 0, 0])
standoff = np.array([0, 20, 0])

# Create coil plane pairs
coil_plus = trimesh.Trimesh(
    planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
)

coil_minus = trimesh.Trimesh(
    planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
)

joined_planes = combine_meshes((coil_plus, coil_minus))

# Create MeshConductor object, which finds the holes and sets the boundary condition
coil = MeshConductor(mesh_obj=joined_planes, fix_normals=True)

#%%
# Set up target and stray field points

# Here, the target points are on a volumetric grid within a sphere

center = np.array([0, 0, 0])

sidelength = 10
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


#%%
# Create bfield specifications used when optimizing the coil geometry

# The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 0] = target_field[:, 0] + 1

target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 0] += 0.001
target_abs_error[:, 1:3] += 0.005

target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": target_abs_error,
    "target": target_field,
}

bfield_specification = [target_spec]

#%%
# Run QP solver
import mosek

coil.s, prob = optimize_streamfunctions(
    coil,
    bfield_specification,
    objective="minimum_inductive_energy",
    solver="MOSEK",
    solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
)

#%%
# Plot the computed streamfunction

coil.s.plot(ncolors=256)
