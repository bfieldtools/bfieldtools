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


font = {"family": "normal", "weight": "normal", "size": 20}
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

Nfreqs = 70
freqs = np.linspace(0, 200, Nfreqs)

#%%
mesh = trimesh.load(
    pkg_resources.resource_filename(
        "bfieldtools", "example_meshes/unitdisc_extremelyfine.stl"
    )
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
Np = 25
z = np.linspace(0.05, 1, Np)
fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T

B_coupling = magnetic_field_coupling(mesh, fp, analytic=True)

Bf = np.sqrt(noise_var(B_coupling, vl))

plt.figure(figsize=(5, 5))
plt.loglog(freqs, Bf[:, 2, :].T * 1e15, linewidth=2)
plt.grid()
plt.ylim(0.1, 100)
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False)
plt.xlabel("Frequency (Hz)")
plt.ylabel(r"$B_z$ noise (fT/rHz)")
plt.tight_layout()

f_interp = np.linspace(0, 120, 500)

cutf = np.zeros(Np)
for i in range(Np):
    Btemp = np.interp(f_interp, freqs, Bf[i, 2, :])
    idx = np.max(np.where(Btemp >= 1 / np.sqrt(2) * Btemp[0]))
    cutf[i] = f_interp[idx]

cutf_an = 1 / (4 * mu0 * sigma * d * z)

plt.figure(figsize=(5, 5))
plt.loglog(z, cutf_an, linewidth=2, label="Infinite plane")
plt.loglog(z, cutf, "x", markersize=10, markeredgewidth=2, label="Disc")
plt.grid()
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.legend(frameon=False, loc=3)
plt.xlabel("Distance (z/R)")
plt.ylabel("3-dB cutoff frequency (Hz)")
plt.tight_layout()
