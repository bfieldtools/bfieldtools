"""
Biplanar coil design
====================

First example in the paper, showing a basic biplanar coil producing homogeneous field in a target
region between the two coil planes.

"""

import pyvista as pv
from bfieldtools._example_plots import point_cloud, set_view

pv.set_plot_theme("document")

PLOT = True
SAVE_FIGURES = False

SAVE_DIR = "./Biplanar coil/"


import numpy as np
import trimesh

from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops
from bfieldtools.utils import load_example_mesh, combine_meshes


# Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


# Load simple plane mesh that is centered on the origin
planemesh = load_example_mesh("10x10_plane_hires")


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

joined_planes = combine_meshes((coil_plus, coil_minus))

joined_planes = joined_planes.subdivide()

# Create mesh class object
coil = MeshConductor(
    verts=joined_planes.vertices,
    tris=joined_planes.faces,
    fix_normals=True,
    basis_name="suh",
    N_suh=100,
)

#%%
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

stray_points_mesh = trimesh.creation.icosphere(subdivisions=3, radius=stray_radius)
stray_points = stray_points_mesh.vertices + center

n_stray_points = len(stray_points)


#%%
# Plot geometry
if PLOT:
    plotter = f = pv.Plotter(window_size=(800, 800))
    plotter = coil.plot_mesh(style="wireframe", opacity=0.1, color=(0, 0, 0), figure=f)
    plotter = coil.plot_mesh(style="surface", opacity=0.1, color=(0, 0, 0), figure=f)
    plotter.add_mesh(
        point_cloud(*target_points.T).glyph(
            scale=False, geom=pv.Sphere(radius=0.05 / 2)
        ),
        color=(0, 0, 1),
        show_scalar_bar=False,
    )
    plotter.add_mesh(
        point_cloud(*stray_points.T).glyph(scale=False, geom=pv.Sphere(radius=0.3 / 2)),
        color=(1, 0, 0),
        show_scalar_bar=False,
    )
    f.view_isometric()
    f.camera.zoom(1.5)
    if SAVE_FIGURES:
        f.screenshot(SAVE_DIR + "biplanar_geometry.png", scale=4)
        plotter.close()


#%%
# Create bfield specifications used when optimizing the coil geometry

# The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 0] += 1  # Homogeneous field on X-axis


target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": 0.005,
    "target": target_field,
}
stray_spec = {
    "coupling": coil.B_coupling(stray_points),
    "abs_error": 0.01,
    "target": np.zeros((n_stray_points, 3)),
}

bfield_specification = [target_spec, stray_spec]

#%%
# Run QP solver

coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec, stray_spec],
    objective=(0, 1),
    solver="CLARABEL",
)

#%%
# Plot coil windings and target points


N_contours = 6

loops = scalar_contour(coil.mesh, coil.s.vert, N_contours=N_contours)

if PLOT:
    plotter = f = pv.Plotter(window_size=(650, 750))
    plotter.clear()
    plotter = plot_3d_current_loops(loops, colors="auto", figure=f)
    f.view_isometric()
    if SAVE_FIGURES:
        f.screenshot(SAVE_DIR + "biplanar_loops.png", scale=4)
        plotter.close()


#%%
# Plot continuous stream function

if PLOT:
    from bfieldtools.viz import plot_data_on_vertices

    plotter = f = pv.Plotter(window_size=(800, 800))
    plotter.clear()
    plotter = plot_data_on_vertices(coil.mesh, coil.s.vert, figure=f, ncolors=256)
    f.camera.parallel_projection = 1
    set_view(plotter, 90, 90)
    f.camera.zoom(1.5)
    if SAVE_FIGURES:
        f.screenshot(SAVE_DIR + "biplanar_streamfunction.png", scale=4)
        plotter.close()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
