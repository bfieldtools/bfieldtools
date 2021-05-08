import numpy as np
import matplotlib.pyplot as plt
import trimesh
from mayavi import mlab
from scipy.linalg import eigh

from bfieldtools.thermal_noise import (
    compute_current_modes,
    visualize_current_modes,
    noise_covar,
    noise_var,
    sensornoise_covar,
)
from bfieldtools.mesh_magnetics import magnetic_field_coupling


import pkg_resources

import mne


#%%

font = {"family": "normal", "weight": "normal", "size": 16}

plt.rc("font", **font)

# Fix the simulation parameters
d = 5e-3
sigma = 3.8e7
T = 293
kB = 1.38064852e-23
mu0 = 4 * np.pi * 1e-7


Nchunks = 4
quad_degree = 2


#%% Save directory of the MNE example here
SAVE_DIR = ""
#%%

evoked = mne.Evoked(SAVE_DIR + "left_auditory-ave.fif")
evoked.pick_types(meg="mag")

coils = mne.forward._create_meg_coils(evoked.info["chs"], acc=2)

#%% Sensor locations and directions in DEVICE coordinate system

p = np.zeros((len(coils), coils[0]["rmag"].shape[0], 3))
n = np.zeros((len(coils), coils[0]["rmag"].shape[0], 3))
w = np.zeros((len(coils), coils[0]["rmag"].shape[0]))

for i in range(len(coils)):
    p[i] = coils[i]["rmag"]
    n[i] = coils[i]["cosmag"]
    w[i] = coils[i]["w"]

#%%

mesh = trimesh.load(
    pkg_resources.resource_filename(
        "bfieldtools", "example_meshes/closed_cylinder_remeshed.stl"
    )
)


theta = np.pi / 2

R = np.zeros((4, 4))
R[0] = np.array((np.cos(theta), 0, np.sin(theta), 0))
R[1] = np.array((0, 1, 0, 0))
R[2] = np.array((-1 * np.sin(theta), 0, np.cos(theta), 0))
R[3] = np.array((0, 0, 0, 0))

mesh.apply_transform(R)

mesh.apply_translation((0, 0, -60e-2))


scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(0, 0, 0), opacity=0.2)
mlab.quiver3d(
    *p.reshape(p.shape[0] * p.shape[1], p.shape[2]).T,
    *n.reshape(n.shape[0] * n.shape[1], n.shape[2]).T,
    mode="arrow"
)
mlab.view(azimuth=0, elevation=0)

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(0, 0, 0), opacity=0.2)

mlab.quiver3d(*p.T, *n.T, mode="arrow")
mlab.view(azimuth=0, elevation=90, distance=5)

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, color=(0, 0, 0), opacity=0.2)

mlab.quiver3d(*p.T, *n.T, mode="arrow")
mlab.view(azimuth=90, elevation=90, distance=5)
#%%

vl, u = compute_current_modes(
    obj=mesh, T=T, resistivity=1 / sigma, thickness=d, mode="AC", return_eigenvals=True
)

noisecov = sensornoise_covar(mesh, p, n, w, vl)
noisecov = noisecov[:, :, 0]

noisecov += noisecov.T - np.diag(np.diag(noisecov))

#%%
plt.figure(figsize=(7.5, 7.5))
plt.imshow(noisecov * 1e30, cmap="RdBu_r")
cbar = plt.colorbar()
plt.clim((-1100, 1100))
plt.xlabel("Sensor index")
plt.ylabel("Sensor index")
cbar.set_label(r"Noise covariance (fT$^2$/Hz)")
plt.tight_layout()

plt.figure(figsize=(2.5, 2.5))
evoked1 = evoked.copy()
evoked1.data = np.tile(np.sqrt(np.diag(noisecov)), (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag", vmin=28, vmax=33)
