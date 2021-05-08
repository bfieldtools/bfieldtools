import numpy as np
import matplotlib.pyplot as plt
import trimesh
from mayavi import mlab
from bfieldtools.thermal_noise import (
    compute_current_modes,
    noise_covar,
    noise_covar_dir,
    noise_var,
    visualize_current_modes,
)
from bfieldtools.mesh_magnetics import magnetic_field_coupling

import pkg_resources


font = {"family": "normal", "weight": "normal", "size": 16}

plt.rc("font", **font)

# Fix the simulation parameters
d = 1e-3
sigma = 3.8e7
T = 293
kB = 1.38064852e-23
mu0 = 4 * np.pi * 1e-7

Nchunks = 4
quad_degree = 2

freqs = np.array((0, 3000))


mesh = trimesh.load(
    pkg_resources.resource_filename("bfieldtools", "example_meshes/star_dense.stl")
)

vl = compute_current_modes(
    obj=mesh,
    T=T,
    resistivity=1 / sigma,
    thickness=d,
    mode="AC",
    freqs=freqs,
    return_eigenvals=False,
)


#%%
scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
Nmodes = 20
dist = 0.2

colormap = "bwr"

N1 = np.floor(np.sqrt(Nmodes))
dx = (mesh.vertices[:, 0].max() - mesh.vertices[:, 0].min()) * (1 + dist)
dy = (mesh.vertices[:, 1].max() - mesh.vertices[:, 1].min()) * (1 + dist)

i = 0
j = 0
for n in range(Nmodes):
    print(i, j)
    points = mesh.vertices.copy()
    points[:, 0] += i * dx
    points[:, 1] -= j * dy

    if n == 0:
        s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(0, 0, 0))
        scene.scene.z_minus_view()
        surface = scene.children[0].children[0].children[0].children[0]
        surface.actor.property.representation = "wireframe"
    else:

        s = mlab.triangular_mesh(
            *points.T,
            mesh.faces,
            scalars=vl[:, n - 1, 0],
            colormap=colormap,
            line_width=6
        )

        limit = np.max(np.abs(vl[:, n]))

        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit, limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.enable_contours = True
        s.contour.number_of_contours = 10

        mlab.triangular_mesh(*points.T, mesh.faces, color=(0, 0, 0), opacity=0.3)

    if i < N1:
        i += 1
    else:
        j += 1
        i = 0

# scene = engine.scenes[0]
scene.scene.z_plus_view()
scene.scene.camera.position = [
    2.4003626448678688,
    -1.7934540101222551,
    12.329167780858096,
]
scene.scene.camera.focal_point = [2.4003626448678688, -1.7934540101222551, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [11.548461583553994, 13.341654730403269]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
scene.scene.camera.position = [
    2.4003626448678688,
    -1.7934540101222551,
    10.189394860213302,
]
scene.scene.camera.focal_point = [2.4003626448678688, -1.7934540101222551, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [9.544183126904127, 11.026160934217575]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()


#%%
Ngrid = 21
xx = np.linspace(np.min(mesh.vertices[:, 0]), np.max(mesh.vertices[:, 0]), Ngrid)
yy = np.linspace(np.min(mesh.vertices[:, 1]), np.max(mesh.vertices[:, 1]), Ngrid)
zz = np.max(mesh.vertices[:, 2]) + np.array([0.1, 0.2, 0.3, 0.4, 0.5])
X, Y, Z = np.meshgrid(xx, yy, zz, indexing="ij")

x = X.ravel()
y = Y.ravel()
z = Z.ravel()

fp = np.vstack((x, y, z)).T


B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)

Bvar = noise_var(B_coupling, vl) * 1e30

#%%

cmap = "inferno"


vmin = 1e-2
vmax = 30


levels = np.linspace(vmin, vmax, 5)

Bzcov = np.sqrt(Bvar[:, 0, 0])
Bgrid = Bzcov.reshape((Ngrid, Ngrid, 5))


from matplotlib import colors
import matplotlib.gridspec as gridspec

Bzcov = np.sqrt(Bvar[:, 0, 0])
Bgrid = Bzcov.reshape((Ngrid, Ngrid, 5))

fig = plt.figure(figsize=(7, 2))
gs = gridspec.GridSpec(1, 5)
gs.update(wspace=0.025, hspace=0.025, left=0, right=0.9, top=1, bottom=0)


for ax_idx in range(5):
    ax = plt.subplot(gs[ax_idx])
    cont = ax.pcolormesh(
        X[:, :, ax_idx],
        Y[:, :, ax_idx],
        Bgrid[:, :, ax_idx],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=colors.LogNorm(),
        shading="gouraud",
    )
    ax.axis("equal")
    ax.axis("off")


cbar_ax = fig.add_axes([0.905, 0.2, 0.01, 0.6])
cbar = fig.colorbar(cont, cax=cbar_ax)


Bzcov = np.sqrt(Bvar[:, 1, 0])
Bgrid = Bzcov.reshape((Ngrid, Ngrid, 5))

fig = plt.figure(figsize=(7, 2))
gs = gridspec.GridSpec(1, 5)
gs.update(wspace=0.025, hspace=0.025, left=0, right=0.9, top=1, bottom=0)

for ax_idx in range(5):
    ax = plt.subplot(gs[ax_idx])
    cont = ax.pcolormesh(
        X[:, :, ax_idx],
        Y[:, :, ax_idx],
        Bgrid[:, :, ax_idx],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=colors.LogNorm(),
        shading="gouraud",
    )
    ax.axis("equal")
    ax.axis("off")


cbar_ax = fig.add_axes([0.905, 0.15, 0.01, 0.7])
cbar = fig.colorbar(cont, cax=cbar_ax)


Bzcov = np.sqrt(Bvar[:, 2, 0])
Bgrid = Bzcov.reshape((Ngrid, Ngrid, 5))

fig = plt.figure(figsize=(7, 2))
gs = gridspec.GridSpec(1, 5)
gs.update(wspace=0.025, hspace=0.025, left=0, right=0.9, top=1, bottom=0)

for ax_idx in range(5):
    ax = plt.subplot(gs[ax_idx])
    cont = ax.pcolormesh(
        X[:, :, ax_idx],
        Y[:, :, ax_idx],
        Bgrid[:, :, ax_idx],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=colors.LogNorm(),
        shading="gouraud",
    )
    ax.axis("equal")
    ax.axis("off")


cbar_ax = fig.add_axes([0.905, 0.2, 0.01, 0.6])
cbar = fig.colorbar(cont, cax=cbar_ax)


Bzcov = np.sqrt(Bvar[:, 2, 1])
Bgrid = Bzcov.reshape((Ngrid, Ngrid, 5))

fig = plt.figure(figsize=(7, 2))
gs = gridspec.GridSpec(1, 5)
gs.update(wspace=0.025, hspace=0.025, left=0, right=0.9, top=1, bottom=0)

for ax_idx in range(5):
    ax = plt.subplot(gs[ax_idx])
    cont = ax.pcolormesh(
        X[:, :, ax_idx],
        Y[:, :, ax_idx],
        Bgrid[:, :, ax_idx],
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        norm=colors.LogNorm(),
        shading="gouraud",
    )
    ax.axis("equal")
    ax.axis("off")


cbar_ax = fig.add_axes([0.905, 0.2, 0.01, 0.6])
cbar = fig.colorbar(cont, cax=cbar_ax)
