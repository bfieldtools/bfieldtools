"""
Magnetically shielded  coil
===========================
Compact example of design of a biplanar coil within a cylindrical shield.
The effect of the shield is prospectively taken into account while designing the coil.
The coil is positioned close to the end of the shield to demonstrate the effect
"""


import pyvista as pv
from bfieldtools._example_plots import vector_glyphs, point_cloud

pv.set_plot_theme("document")


PLOT = True
SAVE_FIGURES = False
SAVE_STREAMFUNCTION = False

SAVE_DIR = "./Shielded coil/"

import numpy as np
import trimesh


from bfieldtools.mesh_conductor import MeshConductor, StreamFunction
from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.utils import load_example_mesh, combine_meshes

from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops

from bfieldtools import resources


# Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


# Load simple plane mesh that is centered on the origin
planemesh = load_example_mesh("10x10_plane_hires")
planemesh.apply_scale(scaling_factor)

# Specify coil plane geometry
center_offset = np.array([9, 0, 0]) * scaling_factor
standoff = np.array([0, 4, 0]) * scaling_factor

# Create coil plane pairs
coil_plus = trimesh.Trimesh(
    planemesh.vertices + center_offset + standoff, planemesh.faces, process=False
)

coil_minus = trimesh.Trimesh(
    planemesh.vertices + center_offset - standoff, planemesh.faces, process=False
)

joined_planes = combine_meshes((coil_minus, coil_plus))


# Create mesh class object
coil = MeshConductor(mesh_obj=joined_planes, fix_normals=True, basis_name="inner")

# Separate object for shield geometry
shieldmesh = trimesh.load(
    file_obj=resources.files(
        "bfieldtools", "example_meshes/closed_cylinder_remeshed.stl"
    ),
    process=True,
)
shieldmesh.apply_scale(15)

shield = MeshConductor(
    mesh_obj=shieldmesh, process=True, fix_normals=True, basis_name="vertex"
)


#%%
# Set up target  points and plot geometry

# Here, the target points are on a volumetric grid within a sphere
# Set up target and stray field points

# Here, the target points are on a volumetric grid within a sphere

center = np.array([9, 0, 0]) * scaling_factor

sidelength = 3 * scaling_factor
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


# Plot coil, shield and target points
if PLOT:
    plotter = f = pv.Plotter(window_size=(800, 800))
    plotter = coil.plot_mesh(style="surface", figure=f, opacity=0.5)
    plotter = shield.plot_mesh(style="surface", opacity=0.2, figure=f)
    plotter.add_mesh(
        point_cloud(*target_points.T).glyph(
            scale=False, geom=pv.Sphere(radius=0.05 / 2)
        ),
        show_scalar_bar=False,
    )
    f.view_isometric()
    f.camera.zoom(1.1)
    if SAVE_FIGURES:
        f.screenshot(SAVE_DIR + "shielded_biplanar_geometry.png", scale=4)
        plotter.close()


#%%
# Let's design a coil without taking the magnetic shield into account

# The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function
target_field = np.zeros(target_points.shape)
target_field[:, 0] = target_field[:, 0] + 1  # Homogeneous Y-field


target_abs_error = np.zeros_like(target_field)
target_abs_error[:, 0] += 0.005
target_abs_error[:, 1:3] += 0.01

target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": target_abs_error,
    "target": target_field,
}


coil.s, coil.prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_inductive_energy",
    solver="CLARABEL",
)


if SAVE_STREAMFUNCTION:
    np.save(SAVE_DIR + "biplanar_streamfunction.npy", coil.s.vert)
#%%
# Plot coil windings and target points

loops = scalar_contour(coil.mesh, coil.s.vert, N_contours=10)

if PLOT:
    plotter = f = pv.Plotter(window_size=(800, 800))
    plotter.clear()
    plotter = plot_3d_current_loops(loops, colors="auto", figure=f)
    B_target = coil.B_coupling(target_points) @ coil.s
    plotter.add_mesh(
        vector_glyphs(*target_points.T, *B_target.T, scale_factor=0.75, mode="arrow"),
        show_scalar_bar=False,
    )
    f.view_isometric()
    f.camera.zoom(0.95)
    if SAVE_FIGURES:
        f.screenshot(SAVE_DIR + "shielded_biplanar_ignored.png", scale=4)
        plotter.close()
#%%
# Now, let's compute the effect of the shield on the field produced by the coil

# Points slightly inside the shield
d = (
    np.mean(np.diff(shield.mesh.vertices[shield.mesh.faces[:, 0:2]], axis=1), axis=0)
    / 10
)
points = shield.mesh.vertices - d * shield.mesh.vertex_normals


# Solve equivalent stream function for the perfect linear mu-metal layer.
# This is the equivalent surface current in the shield that would cause its
# scalar magnetic potential to be constant
shield.s = StreamFunction(
    np.linalg.solve(shield.U_coupling(points), coil.U_coupling(points) @ coil.s), shield
)

#%%
# Plot the difference in field when taking the shield into account

if PLOT:
    plotter = f = pv.Plotter(window_size=(800, 800))
    plotter.clear()
    B_target = coil.B_coupling(target_points) @ coil.s
    B_target_w_shield = (
        coil.B_coupling(target_points) @ coil.s
        + shield.B_coupling(target_points) @ shield.s
    )
    B_quiver = plotter.add_mesh(
        vector_glyphs(
            *target_points.T, *(B_target_w_shield - B_target).T, mode="arrow"
        ),
        cmap="viridis",
        show_scalar_bar=False,
    )
    f.view_isometric()
    plotter.add_scalar_bar(title="Difference in magnetic field (a.u.)")

#%%
# Let's redesign the coil taking the shield into account prospectively

shield.coupling = np.linalg.solve(shield.U_coupling(points), coil.U_coupling(points))

secondary_C = shield.B_coupling(target_points) @ shield.coupling

total_C = coil.B_coupling(target_points) + secondary_C

target_spec_w_shield = {
    "coupling": total_C,
    "abs_error": target_abs_error,
    "target": target_field,
}


coil.s2, coil.prob2 = optimize_streamfunctions(
    coil,
    [target_spec_w_shield],
    objective="minimum_inductive_energy",
    solver="CLARABEL",
)

#%%
# Plot the newly designed coil windings and field at the target points

loops = scalar_contour(coil.mesh, coil.s2.vert, N_contours=10)

if PLOT:
    plotter = f = pv.Plotter(window_size=(800, 800))
    plotter.clear()
    plotter = plot_3d_current_loops(loops, colors="auto", figure=f)
    B_target2 = total_C @ coil.s2
    plotter.add_mesh(
        vector_glyphs(*target_points.T, *B_target2.T, scale_factor=0.75, mode="arrow"),
        show_scalar_bar=False,
    )
    f.view_isometric()
    f.camera.zoom(0.95)
    if SAVE_FIGURES:
        f.screenshot(SAVE_DIR + "shielded_biplanar_prospective.png", scale=4)
        plotter.close()

#%%
# Plot difference in field


import seaborn as sns
import matplotlib.pyplot as plt


if PLOT:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3))

    axnames = ["X", "Y", "Z"]

    # fig.suptitle('Component-wise effect of magnetic shield on target field amplitude distribution')
    for ax_idx, ax in enumerate(axes):

        sns.kdeplot(
            B_target[:, ax_idx],
            label="Coil without shield",
            ax=ax,
            shade=True,
            legend=False,
        )
        sns.kdeplot(
            B_target_w_shield[:, ax_idx],
            label="Coil with shield",
            ax=ax,
            shade=True,
            legend=False,
        )
        sns.kdeplot(
            B_target2[:, ax_idx],
            label="Coil designed with shield",
            ax=ax,
            shade=True,
            legend=False,
        )
        #    ax.set_title(axnames[ax_idx])
        ax.get_yaxis().set_visible(False)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        ax.set_xlabel("Magnetic field on %s-axis" % axnames[ax_idx])

        if ax_idx == 0:
            ax.legend()

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    if SAVE_FIGURES:
        plt.savefig(SAVE_DIR + "shielding_effect.pdf")

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
