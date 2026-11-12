"""
Head gradient coil
==================

Example showing a gradient coil designed on the surface of a MEG system helmet
"""

import pyvista as pv
from bfieldtools._example_plots import vector_glyphs

pv.set_plot_theme("document")


import numpy as np

from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.utils import load_example_mesh
from bfieldtools import sphtools


# Load simple plane mesh that is centered on the origin
helmetmesh = load_example_mesh("meg_helmet")
coil = MeshConductor(mesh_obj=helmetmesh, fix_normals=True)

#%%
# Set up target and stray field points.
# Here, the target points are on a volumetric grid within a sphere

offset = np.array([0, 0, 0.04])
center = offset

sidelength = 0.05
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
# and inner points
target_points = (
    target_points[
        (np.linalg.norm(target_points, axis=1) < sidelength / 2)
        * (np.linalg.norm(target_points, axis=1) > sidelength / 2 * 0.8)
    ]
    + center
)


#%%
# Specify target field and run solver.
# Here, we specify the target field through the use of spherical harmonics.
# We want to produce the field corresponding to a specific beta_l,m-component.

lmax = 3
alm = np.zeros((lmax * (lmax + 2),))
blm = np.zeros((lmax * (lmax + 2),))

# Set one specific component to one
blm[3] += 1

sphfield = sphtools.field(target_points, alm, blm, lmax)

target_field = sphfield / np.max(sphfield[:, 0])

target_field[:, 2] = 0

plotter = coil.plot_mesh(opacity=0.5)
plotter.add_mesh(vector_glyphs(*target_points.T, *sphfield.T), show_scalar_bar=False)
plotter.view_isometric()

abs_error = np.zeros_like(target_field)
abs_error[:, 0] += 0.05
abs_error[:, 1:3] += 0.1


target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": abs_error,
    "target": target_field,
}

#%%

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_inductive_energy",
    solver="CLARABEL",
)

#%%
# Plot coil windings


loops = coil.s.discretize(N_contours=10)
loops.plot_loops()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
