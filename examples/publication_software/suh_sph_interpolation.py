"""
SUH-SPH interpolation comparison
==================================
"""


import numpy as np
from bfieldtools.mesh_conductor import MeshConductor, StreamFunction
from mayavi import mlab
import trimesh
import matplotlib.pyplot as plt

from bfieldtools.sphtools import basis_fields as sphfield
from bfieldtools.sphtools import field as sph_field_eval
from bfieldtools.sphtools import basis_potentials, potential
import mne

from bfieldtools.viz import plot_data_on_vertices, plot_mesh

#%%
SAVE_DIR = "./MNE interpolation/"

#%%
EVOKED = True

with np.load(SAVE_DIR + "mne_data.npz", allow_pickle=True) as data:
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

lmax = 7  # maximum degree
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
PINV = True
if PINV:
    alpha = np.linalg.pinv(Bca_sensors, rcond=1e-15) @ field
else:
    # Calculate using regularization
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

c = MeshConductor(mesh_obj=mesh, basis_name="suh", N_suh=150)
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

f = mlab.figure(bgcolor=(1, 1, 1))
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
fig, ax = plt.subplots(1, 1)
ax.plot(alpha ** 2)


L = np.zeros((0,))
M = np.zeros((0,))


for l in range(1, lmax + 1):
    m_l = np.arange(-l, l + 1, step=1, dtype=np.int_)
    M = np.append(M, m_l)
    L = np.append(L, np.repeat(l, len(m_l)))

xticknames = [None] * len(alpha)
for i in range(len(alpha)):
    xticknames[i] = str(M[i])

    m_l = np.arange(-L[i], L[i] + 1, step=1)

    if i == int(np.floor(len(m_l))):
        xticknames[i] += "\n" + str(L[i])


plt.figure()
plt.plot(a ** 2)


#%% Compute potential on the helmet mesh
from bfieldtools.utils import load_example_mesh
from bfieldtools.flatten_mesh import flatten_mesh, mesh2plane

helmet = load_example_mesh("meg_helmet", process=False)
# Bring the surface roughly to the correct place
helmet.vertices[:, 2] -= 0.045
# The helmet is slightly tilted, correct for this
# (probably the right coordinate transformation could be found from MNE)
rotmat = np.eye(3)
tt = 0.015 * np.pi
rotmat[:2, :2] = np.array([[np.cos(tt), np.sin(tt)], [-np.sin(tt), np.cos(tt)]])
helmet.vertices = helmet.vertices @ rotmat
tt = -0.02 * np.pi
rotmat[1:, 1:] = np.array([[np.cos(tt), np.sin(tt)], [-np.sin(tt), np.cos(tt)]])
helmet.vertices = helmet.vertices @ rotmat
helmet.vertices[:, 1] += 0.005

# plot_mesh(helmet)
# mlab.points3d(*p.T, scale_factor=0.01)


B_sph_helmet = sph_field_eval(
    helmet.vertices,
    alpha,
    np.zeros(alpha.shape),
    lmax=lmax,
    normalization="energy",
    R=R,
)
B_sph_helmet = np.einsum("ij,ij->i", B_sph_helmet, helmet.vertex_normals)
B_suh_helmet = c.B_coupling(helmet.vertices) @ s
B_suh_helmet = np.einsum("ij,ij->i", B_suh_helmet, helmet.vertex_normals)

#%% Compute flattened mesh


u, v, helmet2d = flatten_mesh(helmet, 0.9)
puv = mesh2plane(p, helmet, u, v)


#%% Magnetic field at sensor array surface

from scipy.interpolate import Rbf

rbf_f = Rbf(puv[:, 0], puv[:, 1], field, function="linear", smooth=0)
rbf_field = rbf_f(helmet2d.vertices[:, 0], helmet2d.vertices[:, 1])


vmin = -7e-13
vmax = 7e-13
f = plot_data_on_vertices(helmet2d, rbf_field, ncolors=15, vmin=vmin, vmax=vmax)
mlab.points3d(puv[:, 0], puv[:, 1], 0 * puv[:, 0], scale_factor=0.1, color=(0, 0, 0))
f.scene.z_plus_view()
mlab.savefig(SAVE_DIR + "rbf_helmet_B.png", figure=f, magnification=4)

suh_field = (
    np.einsum("ijk,ij->ik", c.B_coupling(helmet.vertices), helmet.vertex_normals) @ s
)


f = plot_data_on_vertices(helmet2d, suh_field, ncolors=15, vmin=vmin, vmax=vmax)
mlab.points3d(puv[:, 0], puv[:, 1], 0 * puv[:, 0], scale_factor=0.1, color=(0, 0, 0))
f.scene.z_plus_view()
mlab.savefig(SAVE_DIR + "suh_helmet_B.png", figure=f, magnification=4)


Bca, Bcb = sphfield(helmet.vertices, lmax, normalization="energy", R=R)

# sph-components at sensors
sph_field = np.einsum("ijk,ij->ik", Bca, helmet.vertex_normals) @ alpha


f = plot_data_on_vertices(helmet2d, sph_field, ncolors=15, vmin=vmin, vmax=vmax)
mlab.points3d(puv[:, 0], puv[:, 1], 0 * puv[:, 0], scale_factor=0.1, color=(0, 0, 0))
f.scene.z_plus_view()
mlab.savefig(SAVE_DIR + "sph_helmet_B.png", figure=f, magnification=4)

#%% MNE interpolates using splines or something
#%% Compute potential
# U_sph = potential(
# p, alpha, np.zeros(alpha.shape), lmax=lmax, normalization="energy", R=R
# )
#
# U_suh = c.U_coupling(p) @ s

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


# Mask inside/outside using solid angle
mask = abs(c.U_coupling.matrix.sum(axis=1)) < 1e-6
f = plot_data_on_vertices(plane, U_suh * mask, ncolors=15)
# plot_mesh(mesh, figure=f)
f = plot_data_on_vertices(plane, U_sph * mask, ncolors=15)
# plot_mesh(mesh, figure=f)
f = plot_data_on_vertices(plane, (U_suh - U_sph) * mask, ncolors=15)
plot_mesh(mesh, figure=f)
