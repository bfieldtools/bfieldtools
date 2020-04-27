"""
Spherical harmonics-generating coil design
==========================================

Example showing a basic biplanar coil producing a field profile defined by
spherical harmonics. We use the surface harmonics basis for the stream function,
and optimize the coupling to spherical harmonics components, thus creating a compact
optimization problem that can be solved very quickly.

"""

import numpy as np
import trimesh

from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools.coil_optimize import optimize_streamfunctions
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


# To spice things up, let's distort the planes a bit
joined_planes.vertices = (
    joined_planes.vertices
    - 0.5
    * np.linalg.norm(joined_planes.vertices, axis=1)[:, None]
    * joined_planes.vertex_normals
)

# Create mesh class object
coil = MeshConductor(
    mesh_obj=joined_planes, fix_normals=True, basis_name="suh", N_suh=100
)


##############################################################
# Set up target spherical harmonics components

target_alms = np.zeros((coil.opts["N_sph"] * (coil.opts["N_sph"] + 2),))
target_blms = np.zeros((coil.opts["N_sph"] * (coil.opts["N_sph"] + 2),))

target_blms[3] += 1


##############################################################
# Create bfield specifications used when optimizing the coil geometry


target_spec = {
    "coupling": coil.sph_couplings[1],
    "abs_error": 0.01,
    "target": target_blms,
}


##############################################################
# Run QP solver
import mosek

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_ohmic_power",
    solver="MOSEK",
    solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
)


#############################################################
# Plot coil windings

coil.s.discretize(N_contours=8).plot_loops()
