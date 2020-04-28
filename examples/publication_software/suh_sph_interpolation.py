import numpy as np
from bfieldtools.mesh_conductor import MeshConductor, StreamFunction
from mayavi import mlab
import trimesh
import matplotlib.pyplot as plt

from bfieldtools.sphtools import basis_fields as sphfield
from bfieldtools.sphtools import basis_potentials, potential
import mne

#%%
# from mne.datasets import sample

data_path = "/Users/joonas/Codes/MEGNord2019/example_data/MNE-sample-data"

import os

sample_data_raw_file = os.path.join(
    data_path, "MEG", "sample", "sample_audvis_filt-0-40_raw.fif"
)
raw = mne.io.read_raw_fif(sample_data_raw_file)


events = mne.find_events(raw, stim_channel="STI 014")

event_dict = {
    "auditory/left": 1,
    "auditory/right": 2,
    "visual/left": 3,
    "visual/right": 4,
    "smiley": 5,
    "buttonpress": 32,
}

reject_criteria = dict(
    mag=4000e-15,  # 4000 fT
    grad=4000e-13,  # 4000 fT/cm
    eeg=150e-6,  # 150 µV
    eog=250e-6,
)

epochs = mne.Epochs(
    raw, events, event_id=1, tmin=-0.2, tmax=0.5, reject=reject_criteria, preload=True
)
epochs.pick_types(meg="mag")
evoked = epochs.average()

#%%


i0, i1 = evoked.time_as_index(0.08)[0], evoked.time_as_index(0.09)[0]
field = evoked.data[:, i0:i1].mean(axis=1)

# Read BEM for surface geometry and transform to correct coordinate system
import os.path as op

subject = "sample"
subjects_dir = op.join(data_path, "subjects")
bem_fname = op.join(
    subjects_dir, subject, "bem", subject + "-5120-5120-5120-bem-sol.fif"
)
bem = mne.read_bem_solution(bem_fname)

# Head mesh 0
# Innerskull mesh 2
surf_index = 2

trans_fname = op.join(data_path, "MEG", "sample", "sample_audvis_raw-trans.fif")
trans0 = mne.read_trans(trans_fname)
R = trans0["trans"][:3, :3]
t = trans0["trans"][:3, 3]
# Surface from MRI to HEAD
rr = (bem["surfs"][surf_index]["rr"] - t) @ R
# Surface from HEAD to DEVICE
trans1 = evoked.info["dev_head_t"]
R = trans1["trans"][:3, :3]
t = trans1["trans"][:3, 3]
rr = (rr - t) @ R

mesh = trimesh.Trimesh(rr, bem["surfs"][surf_index]["tris"])
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)

surf_index = 0

R = trans0["trans"][:3, :3]
t = trans0["trans"][:3, 3]
# Surface from MRI to HEAD
rr = (bem["surfs"][surf_index]["rr"] - t) @ R
# Surface from HEAD to DEVICE
R = trans1["trans"][:3, :3]
t = trans1["trans"][:3, 3]
rr = (rr - t) @ R
head = trimesh.Trimesh(rr, bem["surfs"][surf_index]["tris"])

mlab.triangular_mesh(*head.vertices.T, head.faces, color=(0.5, 0.5, 0.5), opacity=0.5)

mesh = head


# Sensor locations and directions in DEVICE coordinate system
p = np.array(
    [
        ch["loc"][:3]
        for ch in evoked.info["chs"]
        if ch["ch_name"][-1] == "1" and ch["ch_name"][:3] == "MEG"
    ]
)
n = np.array(
    [
        ch["loc"][-3:]
        for ch in evoked.info["chs"]
        if ch["ch_name"][-1] == "1" and ch["ch_name"][:3] == "MEG"
    ]
)
# Plot sensor locations and directions
mlab.quiver3d(*p.T, *n.T, mode="arrow")


#%% radius for inner/outer sph

R = np.min(np.linalg.norm(p, axis=1)) - 0.02

#%%

lmax = 6  # maximum degree
Bca, Bcb = sphfield(p, lmax, normalization="energy", R=R)

# sph-components at sensors
Bca_sensors = np.einsum("ijk,ij->ik", Bca, n)
Bcb_sensors = np.einsum("ijk,ij->ik", Bcb, n)


#%% Visualize sph components at the helmet
idx = 20

evoked1 = evoked.copy()
evoked1.data[:, :] = np.tile(Bca_sensors[:, idx].T, (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag", colorbar=False)

evoked1 = evoked.copy()
evoked1.data[:, :] = np.tile(Bcb_sensors[:, idx].T, (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag", colorbar=False)

#%% calculate inner sph-coeffients with pinv
alpha = np.linalg.pinv(Bca_sensors, rcond=1e-15) @ field


#%% reconstruct field in helmet

reco_sph = np.zeros(field.shape)

i = 0
for l in range(1, lmax + 1):
    for m in range(-1 * l, l + 1):
        reco_sph += alpha[i] * Bca_sensors[:, i]
        i += 1


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

mlab.figure()
surf = s.plot(False)
surf.actor.mapper.interpolate_scalars_before_mapping = True
surf.module_manager.scalar_lut_manager.number_of_colors = 16

#%% Plot the evoked and the reconsctructions
evoked1 = evoked.copy()
evoked1.data[:, :] = np.tile(field.T, (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag")

evoked1 = evoked.copy()
evoked1.data[:, :] = np.tile(reco_sph.T, (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag")


evoked1 = evoked.copy()
evoked1.data[:, :] = np.tile(reco_suh.T, (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag")


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

#%%
evoked1 = evoked.copy()
evoked1.data[:, :] = np.tile(U_sph.T, (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag")

evoked1 = evoked.copy()
evoked1.data[:, :] = np.tile(U_suh.T, (evoked.times.shape[0], 1)).T
evoked1.plot_topomap(times=0.080, ch_type="mag")
