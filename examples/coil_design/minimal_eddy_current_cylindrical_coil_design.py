"""
Coil with minimal eddy currents
===============================
Compact example of design of a cylindrical coil surrounded by a RF shield, i.e. a conductive surface.
The effects of eddy currents due to inductive interaction with the shield is minimized
"""

import pyvista as pv
from bfieldtools._example_plots import vector_glyphs, point_cloud

pv.set_plot_theme("document")


import numpy as np
import trimesh


from bfieldtools.mesh_conductor import MeshConductor

from bfieldtools.coil_optimize import optimize_streamfunctions
from bfieldtools.contour import scalar_contour
from bfieldtools.viz import plot_3d_current_loops, plot_data_on_vertices

from bfieldtools import resources


# Set unit, e.g. meter or millimeter.
# This doesn't matter, the problem is scale-invariant
scaling_factor = 1


# Load example coil mesh that is centered on the origin
coilmesh = trimesh.load(
    file_obj=resources.files("bfieldtools", "example_meshes/open_cylinder.stl"),
    process=True,
)

angle = np.pi / 2
rotation_matrix = np.array(
    [
        [np.cos(angle), 0, np.sin(angle), 0],
        [0, 1, 0, 0],
        [-np.sin(angle), 0, np.cos(angle), 0],
        [0, 0, 0, 1],
    ]
)

coilmesh.apply_transform(rotation_matrix)

coilmesh1 = coilmesh.copy()

coilmesh2 = coilmesh.copy()


# Create mesh class object
coil = MeshConductor(
    verts=coilmesh1.vertices * 0.75,
    tris=coilmesh1.faces,
    fix_normals=True,
    basis_name="suh",
    N_suh=400,
)


def alu_sigma(T):
    ref_T = 293  # K
    ref_rho = 2.82e-8  # ohm*meter
    alpha = 0.0039  # 1/K

    rho = alpha * (T - ref_T) * ref_rho + ref_rho

    return 1 / rho


resistivity = 1 / alu_sigma(T=293)  # room-temp Aluminium
thickness = 0.5e-3  # 0.5 mm thick


# Separate object for shield geometry
shield = MeshConductor(
    verts=coilmesh2.vertices.copy() * 1.1,
    tris=coilmesh2.faces.copy(),
    fix_normals=True,
    basis_name="inner",
    resistivity=resistivity,
    thickness=thickness,
)


#%%
# Set up target  points and plot geometry


center = np.array([0, 0, 0])

sidelength = 0.25 * scaling_factor
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
plotter = f = pv.Plotter(window_size=(800, 800))
plotter = coil.plot_mesh(figure=f, opacity=0.2)
plotter = shield.plot_mesh(figure=f, opacity=0.2)
plotter.add_mesh(
    point_cloud(*target_points.T).glyph(scale=False, geom=pv.Sphere(radius=0.05 / 2)),
    show_scalar_bar=False,
)


#%%
# Compute eddy-current coupling

mutual_inductance = coil.mutual_inductance(shield)

# Take into account the field produced by currents induced into the shield
# NB! This expression is for instantaneous step-function switching of coil current, see Eq. 18 in G.N. Peeren, 2003.

shield.M_coupling = np.linalg.solve(-shield.inductance, mutual_inductance.T)
secondary_C = shield.B_coupling(target_points) @ -shield.M_coupling

#%%
# Create bfield specifications used when optimizing the coil geometry

# The absolute target field amplitude is not of importance,
# and it is scaled to match the C matrix in the optimization function

target_field = np.zeros(target_points.shape)
target_field[:, 1] = target_field[:, 1] + 1


target_spec = {
    "coupling": coil.B_coupling(target_points),
    "abs_error": 0.01,
    "target": target_field,
}


from scipy.linalg import eigh

l, U = eigh(shield.resistance, shield.inductance, subset_by_index=(0, 500))


time = [0.001, 0.003, 0.005]
eddy_error = [0.05, 0.01, 0.0025]
# time_decay = U @ np.exp(-l[None, :]*time[:, None]) @ np.pinv(U)

time_decay = np.zeros(
    (len(time), shield.inductance.shape[0], shield.inductance.shape[1])
)

induction_spec = []


Uinv = np.linalg.pinv(U)
for idx, t in enumerate(time):
    time_decay = U @ np.diag(np.exp(-l * t)) @ Uinv
    eddy_coupling = shield.B_coupling(target_points) @ time_decay @ shield.M_coupling
    induction_spec.append(
        {
            "coupling": eddy_coupling,
            "abs_error": eddy_error[idx],
            "rel_error": 0,
            "target": np.zeros_like(target_field),
        }
    )

#%%
# Run QP solver to optimize stream function


coil.s, prob = optimize_streamfunctions(
    coil,
    [target_spec] + induction_spec,
    objective="minimum_inductive_energy",
    solver="CLARABEL",
)

from bfieldtools.mesh_conductor import StreamFunction

shield.induced_s = StreamFunction(shield.M_coupling @ coil.s, shield)

#%%
# Plot coil windings and target points


loops = scalar_contour(coil.mesh, coil.s.vert, N_contours=6)


plotter = f = pv.Plotter(window_size=(600, 500))
plotter.clear()

plotter = plot_3d_current_loops(loops, colors="auto", figure=f, tube_radius=0.005)

B_target = coil.B_coupling(target_points) @ coil.s

plotter.add_mesh(vector_glyphs(*target_points.T, *B_target.T), show_scalar_bar=False)

plotter = shield.plot_mesh(
    style="surface",
    opacity=0.5,
    cull_back=True,
    color=(0.8, 0.8, 0.8),
    figure=f,
)
plotter = shield.plot_mesh(
    style="surface",
    opacity=1,
    cull_front=True,
    color=(0.8, 0.8, 0.8),
    figure=f,
)

f.camera.parallel_projection = 1

f.camera.zoom(1.4)

#%%
# For comparison, let's see how the coils look when we ignore the conducting shield


coil.unshielded_s, coil.unshielded_prob = optimize_streamfunctions(
    coil,
    [target_spec],
    objective="minimum_inductive_energy",
    solver="CLARABEL",
)

shield.unshielded_induced_s = StreamFunction(
    shield.M_coupling @ coil.unshielded_s, shield
)

loops = scalar_contour(coil.mesh, coil.unshielded_s.vert, N_contours=6)

plotter = f = pv.Plotter(window_size=(600, 500))
plotter.clear()

plotter = plot_3d_current_loops(loops, colors="auto", figure=f, tube_radius=0.005)

B_target_unshielded = coil.B_coupling(target_points) @ coil.unshielded_s

plotter.add_mesh(
    vector_glyphs(*target_points.T, *B_target_unshielded.T), show_scalar_bar=False
)

plotter = shield.plot_mesh(
    style="surface",
    opacity=0.5,
    cull_back=True,
    color=(0.8, 0.8, 0.8),
    figure=f,
)
plotter = shield.plot_mesh(
    style="surface",
    opacity=1,
    cull_front=True,
    color=(0.8, 0.8, 0.8),
    figure=f,
)

f.camera.parallel_projection = 1

f.camera.zoom(1.4)


#%%
# Finally, let's compare the time-courses


tmin, tmax = 0, 0.025
Fs = 2000

time = np.linspace(tmin, tmax, int(Fs * (tmax - tmin) + 1))

time_decay = np.zeros(
    (len(time), shield.inductance.shape[0], shield.inductance.shape[1])
)

Uinv = np.linalg.pinv(U)
for idx, t in enumerate(time):
    time_decay[idx] = U @ np.diag(np.exp(-l * t)) @ Uinv


B_t = shield.B_coupling(target_points) @ (time_decay @ shield.induced_s).T

unshieldedB_t = (
    shield.B_coupling(target_points) @ (time_decay @ shield.unshielded_induced_s).T
)

import matplotlib.pyplot as plt


fig, ax = plt.subplots(1, 1, sharex=True, figsize=(8, 4))
ax.plot(
    time * 1e3,
    np.mean(np.linalg.norm(B_t, axis=1), axis=0).T,
    "k-",
    label="Minimized",
    linewidth=1.5,
)
ax.set_ylabel("Transient field amplitude")
ax.semilogy(
    time * 1e3,
    np.mean(np.linalg.norm(unshieldedB_t, axis=1), axis=0).T,
    "k--",
    label="Ignored",
    linewidth=1.5,
)
ax.set_xlabel("Time (ms)")


ax.set_ylim(1e-4, 0.5)
ax.set_xlim(0, 25)


plt.grid(which="both", axis="y", alpha=0.1)

plt.legend()
fig.tight_layout()

ax.vlines([1, 5, 10, 20], 1e-4, 0.5, alpha=0.1, linewidth=3, color="r")

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
