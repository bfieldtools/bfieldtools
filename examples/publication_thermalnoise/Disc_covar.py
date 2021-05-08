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


font = {"family": "normal", "weight": "normal", "size": 18}
plt.rc("font", **font)

# Fix the simulation parameters
d = 1e-3
sigma = 3.8e7
T = 293
kB = 1.38064852e-23
mu0 = 4 * np.pi * 1e-7


# freqs = np.array((0,100, 1000, 10000))
freqs = np.array((0, 50, 100, 500))

Np = 121

z = 0.1
x = np.linspace(-1.25, 1.25, Np)
fp = np.array((x, np.zeros(x.shape), z * np.ones(x.shape))).T

ind_orig = 60

Nchunks = 4
quad_degree = 2

mesh = trimesh.load(
    pkg_resources.resource_filename(
        "bfieldtools", "example_meshes/unitdisc_extremelyfine.stl"
    )
)

B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)

vl = compute_current_modes(
    obj=mesh,
    T=T,
    resistivity=1 / sigma,
    thickness=d,
    mode="AC",
    freqs=freqs,
    return_eigenvals=False,
)


Bcov = noise_covar(B_coupling, vl) * (1e30)

Bdircov = noise_covar_dir(B_coupling, vl) * (1e30)


#%%

plt.figure(figsize=(5, 5))
plt.plot(x, np.diag(Bcov[:, :, 0, 0]), linewidth=2, label="$B_x$")
plt.plot(x, np.diag(Bcov[:, :, 1, 0]), linewidth=2, label="$B_y$")
plt.plot(x, np.diag(Bcov[:, :, 2, 0]), linewidth=2, label="$B_z$")
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False)
plt.xlabel("x (m)")
plt.ylabel("Noise variance (fT$^2$/Hz)")
plt.tight_layout()

plt.figure(figsize=(5, 5))
plt.plot(x, Bcov[ind_orig, :, 0, 0], linewidth=2, label="$B_x$")
plt.plot(x, Bcov[ind_orig, :, 1, 0], linewidth=2, label="$B_y$")
plt.plot(x, Bcov[ind_orig, :, 2, 0], linewidth=2, label="$B_z$")
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False)
plt.xlabel("x (m)")
plt.ylabel("Noise covariance (fT$^2$/Hz)")
plt.tight_layout()

plt.figure(figsize=(5, 5))
plt.plot(x, Bcov[ind_orig, :, 2, 0], linewidth=2, label="DC")
plt.plot(x, Bcov[ind_orig, :, 2, 1], linewidth=2, label="50 Hz")
plt.plot(x, Bcov[ind_orig, :, 2, 2], linewidth=2, label="100 Hz")
plt.plot(x, Bcov[ind_orig, :, 2, 3], linewidth=2, label="500 Hz")
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False)
plt.xlabel("x (m)")
plt.ylabel(r"$B_z$ covariance (fT$^2$/Hz)")
plt.tight_layout()


plt.figure(figsize=(3, 3))
plt.plot(x, Bcov[ind_orig, :, 2, 0] / np.max(Bcov[ind_orig, :, 2, 0]), linewidth=2)
plt.plot(x, Bcov[ind_orig, :, 2, 1] / np.max(Bcov[ind_orig, :, 2, 1]), linewidth=2)
plt.plot(x, Bcov[ind_orig, :, 2, 2] / np.max(Bcov[ind_orig, :, 2, 2]), linewidth=2)
plt.plot(x, Bcov[ind_orig, :, 2, 3] / np.max(Bcov[ind_orig, :, 2, 3]), linewidth=2)
plt.gca().tick_params(labelbottom=False, labelleft=False)
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.axis("off")
plt.tight_layout()


plt.figure(figsize=(5, 5))
plt.plot(x, Bdircov[:, 0, 1, 0], linewidth=2, label=r"$B_x$$-$$B_y$")
plt.plot(x, Bdircov[:, 0, 2, 0], linewidth=2, label=r"$B_x$$-$$B_z$")
plt.plot(x, Bdircov[:, 1, 2, 0], linewidth=2, label=r"$B_y$$-$$B_z$")
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False)
plt.xlabel("x (m)")
plt.ylabel("Noise covariance (fT$^2$/Hz)")
plt.tight_layout()
