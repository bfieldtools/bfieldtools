"""
Thermal noise above a conducting disc
=====================================

Study mesh convergence and current modes for thermal magnetic noise
above a thin conducting disc.
"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh

pv.set_plot_theme("document")

#%%
import numpy as np
import matplotlib.pyplot as plt
import trimesh
from scipy.linalg import eigh

from bfieldtools.thermal_noise import (
    compute_current_modes,
    visualize_current_modes,
    noise_covar,
    noise_var,
)
from bfieldtools.mesh_magnetics import magnetic_field_coupling

from bfieldtools import utils

from bfieldtools import resources

font = {"family": "normal", "weight": "normal", "size": 16}
plt.rc("font", **font)


#%%
# Fix the simulation parameters
d = 1e-3
sigma = 3.8e7
T = 293
kB = 1.38064852e-23
mu0 = 4 * np.pi * 1e-7

freqs = np.array((0,))

Np = 200

z = np.linspace(0.05, 10, Np)
fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T

Niters = 3
Nfaces = np.zeros(Niters)
Bn = np.zeros((Niters, Np))

modinds = np.array((0, 4, 9, 49, 99, 249))
Nmods = len(modinds)

Nchunks = 4
quad_degree = 4

for i in range(Niters):
    if i == 0:
        mesh = trimesh.load(
            resources.files("bfieldtools", "example_meshes/disc_finer.stl")
        )
    if i == 1:
        mesh = trimesh.load(
            resources.files("bfieldtools", "example_meshes/disc_extrafine.stl")
        )
    if i == 2:
        mesh = trimesh.load(
            resources.files("bfieldtools", "example_meshes/disc_extremelyfine.stl")
        )
    B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)
    Nfaces[i] = mesh.faces.shape[0]
    (vl, u) = compute_current_modes(
        obj=mesh,
        T=T,
        resistivity=1 / sigma,
        thickness=d,
        mode="AC",
        return_eigenvals=True,
    )
    plotter = scene = pv.Plotter(window_size=(800, 800))
    visualize_current_modes(mesh, vl[:, :, 0], 20, 5, contours=True)
    for m in range(Nmods):
        plotter = scene = pv.Plotter(window_size=(800, 800))
        s = plotter.add_mesh(
            triangle_mesh(*mesh.vertices.T, mesh.faces),
            scalars=vl[:, modinds[m], 0],
            line_width=12,
            cmap="bwr",
            show_scalar_bar=False,
        )
        limit = np.max(np.abs(vl[:, modinds[m], 0]))
        mesh_actor(s).mapper.scalar_range = np.array([-limit, limit])
        mesh_actor(s).mapper.dataset = mesh_actor(s).mapper.dataset.contour(10)
        scene.view_xy()
        scene.camera.position = [
            -5.513350804725592e-05,
            -1.691800821806977e-05,
            4.020380431659883,
        ]
        scene.camera.focal_point = [-5.513350804725592e-05, -1.691800821806977e-05, 0.0]
        scene.camera.view_angle = 30.0
        scene.camera.up = [0.0, 1.0, 0.0]
        scene.camera.clipping_range = [3.7658023470473827, 4.350539189462033]
        scene.render()
    pv.close_all()
    plotter = scene = pv.Plotter(window_size=(800, 800))
    s = plotter.add_mesh(
        triangle_mesh(*mesh.vertices.T, mesh.faces), show_scalar_bar=False
    )
    scene.view_xy(negative=True)
    s.prop.style = "wireframe"
    s.mapper.scalar_visibility = False
    scene.camera.position = [0.0, 0.0, -4.515786458791532]
    scene.camera.focal_point = [0.0, 0.0, 0.0]
    scene.camera.view_angle = 30.0
    scene.camera.up = [0.0, 1.0, 0.0]
    scene.camera.clipping_range = [4.229838328573525, 4.886628590046963]
    scene.render()
    Bvar = noise_var(B_coupling, vl)
    Bn[i] = np.sqrt(Bvar[:, 2, 0])
    if i == 0:
        u0 = u
    if i == 1:
        u1 = u
    if i == 2:
        u2 = u
    print(i)

figw = 3.75

plt.figure(figsize=(figw, 5))
plt.loglog(1 / u0 * 1e3, linewidth=2, label="N = %i" % Nfaces[0])
plt.loglog(1 / u1 * 1e3, linewidth=2, label="N = %i" % Nfaces[1])
plt.loglog(1 / u2 * 1e3, linewidth=2, label="N = %i" % Nfaces[2])
plt.ylim((0.1, 10))
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.xlabel("Mode index")
plt.ylabel("Time constant (ms)")
plt.legend(frameon=False, loc=3)
plt.tight_layout()


r = 1
Ban = (
    mu0
    * np.sqrt(sigma * d * kB * T / (8 * np.pi * z**2))
    * (1 / (1 + z**2 / r**2))
)

plt.figure(figsize=(figw, 5))
plt.loglog(z, Ban * 1e15, linewidth=2, label="Analytic")
plt.loglog(z, Bn[2] * 1e15, "x", label="N = %i" % Nfaces[2])
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False, loc=3)
plt.xlabel("Distance (z/R)")
plt.ylabel(r"$B_z$ noise at DC (fT/rHz)")
plt.tight_layout()

plt.figure(figsize=(figw, 5))
for i in range(Niters):
    plt.loglog(
        z,
        np.abs((Bn[i] - Ban)) / np.abs(Ban) * 100,
        linewidth=2,
        label="N = %i" % Nfaces[i],
    )
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.xlabel("Distance (z/R)")
plt.ylabel("Relative error (%)")
plt.legend(frameon=False)
plt.tight_layout()

N = 100
Nmodes = np.floor(np.linspace(1, vl.shape[1], N)).astype("int")
Bm = np.zeros((N, Np))

for k in range(N):
    Bvar = noise_var(B_coupling, vl, Nmodes=Nmodes[k])
    Bm[k] = np.sqrt(Bvar[:, 2, 0])


plt.figure(figsize=(figw, 5))
Ban0 = (
    mu0
    * np.sqrt(sigma * d * kB * T / (8 * np.pi * z[0] ** 2))
    * (1 / (1 + z[0] ** 2 / r**2))
)
plt.semilogy(
    Nmodes,
    np.abs((Bm[:, 0] - Ban0)) / np.abs(Ban0) * 100,
    linewidth=2,
    label="z/R = %.2f" % z[0],
)

Ban0 = (
    mu0
    * np.sqrt(sigma * d * kB * T / (8 * np.pi * z[1] ** 2))
    * (1 / (1 + z[1] ** 2 / r**2))
)
plt.semilogy(
    Nmodes,
    np.abs((Bm[:, 1] - Ban0)) / np.abs(Ban0) * 100,
    linewidth=2,
    label="z/R = %.2f" % z[1],
)

plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.xlabel("Number of modes")
plt.ylabel("Relative error (%)")
plt.legend(frameon=False)
plt.tight_layout()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
