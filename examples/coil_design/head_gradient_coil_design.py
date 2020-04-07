"""
Head gradient coil
==================

Example showing a gradient coil designed on the surface of a MEG system helmet
"""


import numpy as np
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
helmetmesh = trimesh.load(
    file_obj=pkg_resources.resource_filename(
        "bfieldtools", "example_meshes/meg_helmet.obj"
    ),
    process=False,
)

# planemesh.apply_scale(scaling_factor)
#
##Specify coil plane geometry
# center_offset = np.array([0, 0, 0]) * scaling_factor
# standoff = np.array([0, 4, 0]) * scaling_factor
#
##Create coil plane pairs
# coil_plus = trimesh.Trimesh(planemesh.vertices + center_offset + standoff,
#                         planemesh.faces, process=False)
##
# coil_minus = trimesh.Trimesh(planemesh.vertices + center_offset - standoff,
#                     planemesh.faces, process=False)

# joined_planes = coil_plus.union(coil_minus)

# Create mesh class object
coil = Conductor(verts=helmetmesh.vertices, tris=helmetmesh.faces, fix_normals=True)

###############################################################
# Set up target and stray field points.
# Here, the target points are on a volumetric grid within a sphere

offset = np.array([0, 0, 0.04])
center = offset * scaling_factor

sidelength = 0.05 * scaling_factor
n = 12
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


###############################################################
# Specify target field and run solver

# Let's generate the target field through the use of spherical harmonics.
# Thus we avoid issues with having to manually specify the concomitant gradients


from bfieldtools import sphtools


lmax = 3
alm = np.zeros((lmax * (lmax + 2),))
blm = np.zeros((lmax * (lmax + 2),))

#

blm[3] += 1

sphfield = sphtools.field(target_points - offset, alm, blm, lmax)

target_field = sphfield / np.max(sphfield[:, 0])

target_field[:, 2] = 0

coil.plot_mesh()
mlab.quiver3d(*target_points.T, *sphfield.T)


rel_error = np.zeros_like(target_field)
# rel_error[:, 0] += 0.1

abs_error = np.zeros_like(target_field)
abs_error[:, 0] += 0.1
abs_error[:, 1:3] += 0.1


target_spec = {
    "coupling": coil.B_coupling(target_points),
    "rel_error": rel_error,
    "abs_error": abs_error,
    "target": target_field,
}

import mosek

coil.j, prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_inductive_energy",
    solver="MOSEK",
    solver_opts={"mosek_params": {mosek.iparam.num_threads: 8}},
)

###############################################################
# Plot coil windings and magnetic field in target points


loops = scalar_contour(coil.mesh, coil.j, N_contours=20)

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
mlab.clf()

plot_3d_current_loops(loops, colors="auto", figure=f, tube_radius=0.05 / 50)

B_target = coil.B_coupling(target_points) @ coil.j

mlab.quiver3d(*target_points.T, *B_target.T)

f.scene.isometric_view()
