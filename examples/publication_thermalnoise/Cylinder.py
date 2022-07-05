import numpy as np
import matplotlib.pyplot as plt
import trimesh
from mayavi import mlab
from bfieldtools.thermal_noise import (
    compute_current_modes,
    noise_covar,
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
freqs = np.array((0,))


Nchunks = 8
quad_degree = 2


mesh = trimesh.load(
    pkg_resources.resource_filename("bfieldtools", "example_meshes/closed_cylinder.stl")
)
mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)


vl = compute_current_modes(
    obj=mesh, T=T, resistivity=1 / sigma, thickness=d, mode="AC", return_eigenvals=False
)


vl[:, 0] = np.zeros(vl[:, 0].shape)  # fix DC-component

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))

visualize_current_modes(mesh, vl[:, :, 0], 8, 1)


scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
scene.scene.z_minus_view()
surface = scene.children[0].children[0].children[0].children[0]
surface.actor.property.representation = "wireframe"
surface.actor.mapper.scalar_visibility = False
scene.scene.isometric_view()
scene.scene.render()


Np = 30

x = np.linspace(-0.95, 0.95, Np)
fp = np.array((x, np.zeros(x.shape), np.zeros(x.shape))).T

B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)
B = noise_var(B_coupling, vl)

a = 0.5
L = 2
rat = L / (2 * a)
Gfact = (
    1
    / (8 * np.pi)
    * (
        (3 * rat**5 + 5 * rat**3 + 2) / (rat**2 * (1 + rat**2) ** 2)
        + 3 * np.arctan(rat)
    )
)
Ban = np.sqrt(Gfact) * mu0 * np.sqrt(kB * T * sigma * d) / a

plt.figure(figsize=(5, 5))
plt.plot(x, Ban * np.ones(x.shape) * 1e15, label="Analytic", linewidth=2)
plt.plot(
    x,
    np.sqrt(B[:, 0]) * 1e15,
    "x",
    label="Numerical",
    markersize=10,
    markeredgewidth=2,
)
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False)
plt.xlabel("z (m)")
plt.ylabel(r"$B_z$ noise at DC (fT/rHz)")
plt.tight_layout()

plt.figure()
plt.semilogy(x, np.sqrt(B[:, 0]), label="x")
plt.semilogy(x, np.sqrt(B[:, 1]), label="y")
plt.semilogy(x, np.sqrt(B[:, 2]), "--", label="z")
plt.legend()
plt.xlabel("z (m)")
plt.ylabel("DC noise (T/rHz)")
