"""
Spherical harmonics-generating coil design
==========================================

Example showing a basic biplanar coil producing a field profile defined by
spherical harmonics.

"""

# import sys
# path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
# if path in sys.path:
#    sys.path.insert(0, path)


import numpy as np
from mayavi import mlab
import trimesh

from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops
from bfieldtools import sphtools


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

planemesh.apply_scale(scaling_factor)

# Specify coil plane geometry
center_offset = np.array([0, 0, 0]) * scaling_factor
standoff = np.array([0, 4, 0]) * scaling_factor

# Create coil plane pairs
coil_plus = trimesh.Trimesh(
    planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
)

coil_minus = trimesh.Trimesh(
    planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
)

joined_planes = coil_plus.union(coil_minus)

# Create mesh class object

coil = MeshConductor(
    verts=joined_planes.vertices,
    tris=joined_planes.faces,
    fix_normals=True,
    N_sph=4,
    basis_name="suh",
    N_suh=80,
)

target_alms = np.zeros((coil.opts["N_sph"] * (coil.opts["N_sph"] + 2),))
target_blms = np.zeros((coil.opts["N_sph"] * (coil.opts["N_sph"] + 2),))

target_blms[3] += 1


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


sphfield = sphtools.field(target_points, target_alms, target_blms, coil.opts["N_sph"])


target_field = sphfield / np.max(np.abs(sphfield))

coil.plot_mesh()
mlab.quiver3d(*target_points.T, *sphfield.T)


##############################################################
# Create bfield specifications used when optimizing the coil geometry


target_spec = {
    "coupling": coil.sph_couplings[1],
    "rel_error": 0,
    "abs_error": 0.01,
    "target": target_blms,
}


##############################################################
# Run QP solver
import mosek

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_inductive_energy",
    solver="MOSEK",
    solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
)


#############################################################
# Plot coil windings and target points

N_contours = 10

loops = scalar_contour(coil.mesh, coil.s.vert, N_contours=N_contours)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors="auto", figure=f)

B_target = coil.B_coupling(target_points) @ coil.s

mlab.quiver3d(*target_points.T, *B_target.T)
