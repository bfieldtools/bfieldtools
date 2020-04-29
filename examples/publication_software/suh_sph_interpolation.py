import numpy as np
from bfieldtools.mesh_conductor import MeshConductor, StreamFunction
from mayavi import mlab
import trimesh
import matplotlib.pyplot as plt

from bfieldtools.sphtools import basis_fields as sphfield
from bfieldtools.sphtools import basis_potentials, potential
import mne


#%%
SAVE_DIR = "./MNE interpolation/"

#%%
EVOKED = True

with np.load(SAVE_DIR + "mne_data.npz", allow_pickle=True) as data:
    mesh = data["mesh"]
    p = data["p"]
    n = data["n"]
    mesh = trimesh.Trimesh(vertices=data["vertices"], faces=data["faces"])

if EVOKED:
    evoked = mne.Evoked(SAVE_DIR + "left_auditory-ave.fif")

    i0, i1 = evoked.time_as_index(0.08)[0], evoked.time_as_index(0.09)[0]
    field = evoked.data[:, i0:i1].mean(axis=1)

else:
    # take "data" from lead field matrix, i.e, topography of a single dipole
    from mne.datasets import sample
    import os

    data_path = sample.data_path()

    raw_fname = data_path + "/MEG/sample/sample_audvis_raw.fif"
    trans = data_path + "/MEG/sample/sample_audvis_raw-trans.fif"
    src = data_path + "/subjects/sample/bem/sample-oct-6-src.fif"
    bem = data_path + "/subjects/sample/bem/sample-5120-5120-5120-bem-sol.fif"
    subjects_dir = os.path.join(data_path, "subjects")

    # Note that forward solutions can also be read with read_forward_solution
    fwd = mne.make_forward_solution(
        raw_fname, trans, src, bem, meg=True, eeg=False, mindist=5.0, n_jobs=2
    )
    # Take only magnetometers
    mags = np.array([n[-1] == "1" for n in fwd["sol"]["row_names"]])
    L = fwd["sol"]["data"][mags, :]
    # Take the first dipole
    field = L[:, 56]

#%% radius for inner/outer sph

R = np.min(np.linalg.norm(p, axis=1)) - 0.02

#%%

lmax = 9  # maximum degree
Bca, Bcb = sphfield(p, lmax, normalization="energy", R=R)

# sph-components at sensors
Bca_sensors = np.einsum("ijk,ij->ik", Bca, n)
Bcb_sensors = np.einsum("ijk,ij->ik", Bcb, n)


#%% Visualize sph components at the helmet
# idx = 20

# evoked1 = evoked.copy()
# evoked1.data[:, :] = np.tile(Bca_sensors[:, idx].T, (evoked.times.shape[0], 1)).T
# evoked1.plot_topomap(times=0.080, ch_type="mag", colorbar=False)

# evoked1 = evoked.copy()
# evoked1.data[:, :] = np.tile(Bcb_sensors[:, idx].T, (evoked.times.shape[0], 1)).T
# evoked1.plot_topomap(times=0.080, ch_type="mag", colorbar=False)

#%% calculate inner sph-coeffients with pinv
PINV = False
if PINV:
    alpha = np.linalg.pinv(Bca_sensors, rcond=1e-15) @ field
else:
    ssa = np.linalg.svd(Bca_sensors @ Bca_sensors.T, False, False)
    reg_exp = 6
    _lambda = np.max(ssa) * (10 ** (-reg_exp))
    # angular-Laplacian in the sph basis is diagonal
    La = np.diag([l * (l + 1) for l in range(1, lmax + 1) for m in range(-l, l + 1)])
    BB = Bca_sensors.T @ Bca_sensors + _lambda * La
    alpha = np.linalg.solve(BB, Bca_sensors.T @ field)

# Reconstruct field in helmet

# reco_sph = np.zeros(field.shape)
# i = 0
# for l in range(1, lmax + 1):
#     for m in range(-1 * l, l + 1):
#         reco_sph += alpha[i] * Bca_sensors[:, i]
#         i += 1

# Produces the same result as the loop
reco_sph = Bca_sensors @ alpha

print(
    "SPH-reconstruction relative error:",
    np.linalg.norm(reco_sph - field) / np.linalg.norm(field),
)

#%%
##%% Fit the surface current for the auditory evoked response using pinv
# c = MeshConductor(mesh_obj=mesh, basis_name="suh", N_suh=35)
# M = c.mass
# B_sensors = np.einsum("ijk,ij->ik", c.B_coupling(p), n)
#
#
# asuh = np.linalg.pinv(B_sensors, rcond=1e-15) @ field
#
# s = StreamFunction(asuh, c)
# b_filt = B_sensors @ s


#%% Suh fit

c = MeshConductor(mesh_obj=mesh, basis_name="suh", N_suh=250)
M = c.mass

B_sensors = np.einsum("ijk,ij->ik", c.B_coupling(p), n)
ss = np.linalg.svd(B_sensors @ B_sensors.T, False, False)

reg_exp = 1
plot_this = True
rel_errors = []
_lambda = np.max(ss) * (10 ** (-reg_exp))
# Laplacian in the suh basis is diagonal
BB = B_sensors.T @ B_sensors + _lambda * (-c.laplacian) / np.max(abs(c.laplacian))
a = np.linalg.solve(BB, B_sensors.T @ field)

s = StreamFunction(a, c)

reco_suh = B_sensors @ s

print(
    "SUH-reconstruction relative error:",
    np.linalg.norm(reco_suh - field) / np.linalg.norm(field),
)

f = mlab.figure()
surf = s.plot(False, figure=f)
surf.actor.mapper.interpolate_scalars_before_mapping = True
surf.module_manager.scalar_lut_manager.number_of_colors = 16

#%% Plot the evoked and the reconsctructions
# evoked1 = evoked.copy()
# evoked1.data[:, :] = np.tile(field.T, (evoked.times.shape[0], 1)).T
# evoked1.plot_topomap(times=0.080, ch_type="mag")

# evoked1 = evoked.copy()
# evoked1.data[:, :] = np.tile(reco_sph.T, (evoked.times.shape[0], 1)).T
# evoked1.plot_topomap(times=0.080, ch_type="mag")


# evoked1 = evoked.copy()
# evoked1.data[:, :] = np.tile(reco_suh.T, (evoked.times.shape[0], 1)).T
# evoked1.plot_topomap(times=0.080, ch_type="mag")


#%% Plot spectra
plt.figure()
plt.plot(alpha ** 2)

plt.figure()
plt.plot(a ** 2)


#%% Compute potential
U_sph = potential(
    p, alpha, np.zeros(alpha.shape), lmax=lmax, normalization="energy", R=R
)

U_suh = c.U_coupling(p) @ a

#%% MNE interpolates using SPHERICAL HEAD MODEL
# evoked1 = evoked.copy()
# evoked1.data[:, :] = np.tile(U_sph.T, (evoked.times.shape[0], 1)).T
# evoked1.plot_topomap(times=0.080, ch_type="mag")

# evoked1 = evoked.copy()
# evoked1.data[:, :] = np.tile(U_suh.T, (evoked.times.shape[0], 1)).T
# evoked1.plot_topomap(times=0.080, ch_type="mag")


#%% interpolate data on planar mesh
from bfieldtools.utils import load_example_mesh
from bfieldtools.mesh_calculus import gradient

plane = load_example_mesh("10x10_plane_hires")
scaling_factor = 0.03
plane.apply_scale(scaling_factor)
# Rotate to x-plane
t = np.eye(4)
theta = np.pi / 2 * 1.2
t[1:3, 1:3] = np.array(
    [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]
)
plane.apply_transform(t)

c.U_coupling.reset()
U_suh = c.U_coupling(plane.vertices) @ a
# Adapt mesh to the function and calculate new points
for i in range(2):
    g = np.linalg.norm(gradient(U_suh, plane), axis=0)
    face_ind = np.flatnonzero(g > g.max() * 0.05)
    plane = plane.subdivide(face_ind)
    U_suh = c.U_coupling(plane.vertices) @ a

U_sph = potential(
    plane.vertices, alpha, np.zeros(alpha.shape), lmax=lmax, normalization="energy", R=R
)

#%%
from bfieldtools.viz import plot_data_on_vertices
from bfieldtools.viz import plot_mesh

# Mask inside/outside using solid angle
mask = abs(c.U_coupling.matrix.sum(axis=1)) < 1e-6
f = plot_data_on_vertices(plane, U_suh * mask, ncolors=15)
# plot_mesh(mesh, figure=f)
f = plot_data_on_vertices(plane, U_sph * mask, ncolors=15)
# plot_mesh(mesh, figure=f)
f = plot_data_on_vertices(plane, (U_suh - U_sph) * mask, ncolors=15)
plot_mesh(mesh, figure=f)
