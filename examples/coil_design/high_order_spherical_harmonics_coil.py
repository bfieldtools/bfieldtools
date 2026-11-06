"""
High-order spherical harmonic biplanar coil design
==================================================

Example showing a basic biplanar coil producing a high-order spherical harmonic field
in a specific target region between the two coil planes.

"""

import pyvista as pv
from bfieldtools._example_plots import vector_glyphs

pv.set_plot_theme("document")


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


#%%
# Create bfield specifications used when optimizing the coil geometry


from bfieldtools import sphtools


lmax = 4
alm = np.zeros((lmax * (lmax + 2),))
blm = np.zeros((lmax * (lmax + 2),))

# Set one specific component to one
blm[16] += 1

sphfield = sphtools.field(target_points, alm, blm, lmax)

target_field = sphfield / np.max(sphfield[:, 0])


plotter = coil.plot_mesh(opacity=0.2)
plotter.add_mesh(vector_glyphs(*target_points.T, *sphfield.T), show_scalar_bar=False)


target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": 0.1,
    "target": target_field,
}


#%%
# Run QP solver

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_inductive_energy",
    solver="CLARABEL",
)

#%%
# Plot coil windings and target points
coil.s.discretize(N_contours=10).plot_loops()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
