"""
Interference rejection example using equivalent surface currents
==================================================================

"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh, vector_glyphs, structured_field

pv.set_plot_theme("document")

import numpy as np
from bfieldtools.mesh_conductor import MeshConductor, StreamFunction
import trimesh

import mne

PLOT = True
SAVE_FIGURES = False
IMPORT_MNE_DATA = True

SAVE_MNE_DATA = True
SAVE_DIR = "./MNE interpolation/"


#%%
# First, let's import the MEG data

if IMPORT_MNE_DATA:
    from mne.datasets import sample

    data_path = sample.data_path()
    fname = data_path / "MEG/sample/sample_audvis-ave.fif"
    condition = "Left Auditory"
    evoked = mne.read_evokeds(fname, condition=condition, baseline=(None, 0), proj=True)
    evoked.pick_types(meg="mag")
    (i0, i1) = (evoked.time_as_index(0.08)[0], evoked.time_as_index(0.09)[0])
    field = evoked.data[:, i0:i1].mean(axis=1)
    import os.path as op

    subject = "sample"
    subjects_dir = op.join(data_path, "subjects")
    bem_fname = op.join(
        subjects_dir, subject, "bem", subject + "-5120-5120-5120-bem-sol.fif"
    )
    bem = mne.read_bem_solution(bem_fname)
    surf_index = 2
    trans_fname = op.join(data_path, "MEG", "sample", "sample_audvis_raw-trans.fif")
    trans0 = mne.read_trans(trans_fname)
    R = trans0["trans"][:3, :3]
    t = trans0["trans"][:3, 3]
    rr = (bem["surfs"][surf_index]["rr"] - t) @ R
    trans1 = evoked.info["dev_head_t"]
    R = trans1["trans"][:3, :3]
    t = trans1["trans"][:3, 3]
    rr = (rr - t) @ R
    innerskull = trimesh.Trimesh(rr, bem["surfs"][surf_index]["tris"])
    surf_index = 0
    R = trans0["trans"][:3, :3]
    t = trans0["trans"][:3, 3]
    rr = (bem["surfs"][surf_index]["rr"] - t) @ R
    R = trans1["trans"][:3, :3]
    t = trans1["trans"][:3, 3]
    rr = (rr - t) @ R
    head = trimesh.Trimesh(rr, bem["surfs"][surf_index]["tris"])
    mesh = head
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
    if PLOT:
        plotter = fig = pv.Plotter()
        plotter.add_mesh(
            triangle_mesh(*innerskull.vertices.T, innerskull.faces),
            show_scalar_bar=False,
        )
        plotter.add_mesh(
            triangle_mesh(*head.vertices.T, head.faces),
            color=(0.5, 0.5, 0.5),
            opacity=0.5,
            show_scalar_bar=False,
        )
        plotter.add_mesh(vector_glyphs(*p.T, *n.T, mode="arrow"), show_scalar_bar=False)
        fig.view_isometric()
        if SAVE_FIGURES:
            fig.screenshot(SAVE_DIR + "MEG_geometry.png", scale=4)
    if SAVE_MNE_DATA:
        np.savez(
            SAVE_DIR + "mne_data.npz",
            mesh=head,
            p=p,
            n=n,
            vertices=head.vertices,
            faces=head.faces,
        )
        evoked.save(SAVE_DIR + "left_auditory-ave.fif")
else:
    with np.load(SAVE_DIR + "mne_data.npz", allow_pickle=True) as data:
        mesh = data["mesh"]
        p = data["p"]
        n = data["n"]
        mesh = trimesh.Trimesh(vertices=data["vertices"], faces=data["faces"])
    evoked = mne.Evoked(SAVE_DIR + "left_auditory-ave.fif")


#%% Construct meshes, compute couplings
# Fit the surface current for the auditory evoked response

SCALING = 3

from bfieldtools.utils import combine_meshes

# mesh2 = trimesh.creation.icosphere(subdivisions=5, radius=2, color=None)

mesh2 = mesh.copy()
mesh2.vertices *= SCALING

cmesh = combine_meshes((mesh, mesh2))

c = MeshConductor(mesh_obj=cmesh)  # , basis_name="suh", N_suh=300)
M = c.mass


sensor_coupling = np.einsum("ijk,ij->ik", c.B_coupling(p), n)
# a = np.linalg.pinv(sensor_coupling, rcond=1e-15) @ field
ss = np.linalg.svd(sensor_coupling @ sensor_coupling.T, False, False)


#%% Fit surface currents without any outside interference

# reg_exps = [0.5, 1, 2, 3, 4, 5, 6, 7, 8]
reg_exps = [1]
rel_errors = []

for reg_exp in reg_exps:
    _lambda = np.max(ss) * 10 ** (-reg_exp)
    BB = sensor_coupling.T @ sensor_coupling + _lambda * -c.laplacian / np.max(
        abs(c.laplacian)
    )
    a = np.linalg.solve(BB, sensor_coupling.T @ field)
    s = StreamFunction(a, c)
    b_filt = sensor_coupling @ s
    rel_error = np.linalg.norm(b_filt - field) / np.linalg.norm(field)
    print("Fit relative error:", rel_error * 100, "%")
    rel_errors.append(rel_error)
    if PLOT:
        plotter = fig = pv.Plotter()
        surf1 = plotter.add_mesh(
            triangle_mesh(*mesh.vertices.T, mesh.faces),
            scalars=s.vert[0 : len(mesh.vertices)],
            show_scalar_bar=False,
        )
        mesh_actor(surf1).mapper.interpolate_before_map = True
        mesh_actor(surf1).mapper.lookup_table.n_values = 16
        surf2 = plotter.add_mesh(
            triangle_mesh(*mesh2.vertices.T, mesh2.faces),
            scalars=s.vert[len(mesh.vertices) :],
            show_scalar_bar=False,
        )
        mesh_actor(surf2).mapper.dataset = mesh_actor(surf2).mapper.dataset.contour(10)
        if SAVE_FIGURES:
            fig.screenshot(SAVE_DIR + "SUH_scalp_streamfunction.png", scale=4)


#%% Add outside interference to data and test fit

from bfieldtools.line_conductor import LineConductor


currentloop = LineConductor(
    [np.array([[10, 1, 1], [10, 1, -1], [10, -1, -1], [10, -1, 1], [10, 1, 1]])]
)


CURRENT = 1e-4


interference = np.einsum("ij,ij->i", (currentloop.magnetic_field(p) * CURRENT), n)

dirtyfield = field + interference


reg_exps = [-1, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]
# reg_exps = [1]
rel_errors = []

for reg_exp in reg_exps:
    _lambda = np.max(ss) * 10 ** (-reg_exp)
    BB = sensor_coupling.T @ sensor_coupling + _lambda * -c.laplacian / np.max(
        abs(c.laplacian)
    )
    a = np.linalg.solve(BB, sensor_coupling.T @ dirtyfield)
    s = StreamFunction(a, c)
    b_filt = sensor_coupling @ s
    rel_error = np.linalg.norm(b_filt - dirtyfield) / np.linalg.norm(dirtyfield)
    print("Fit relative error:", rel_error * 100, "%")
    rel_errors.append(rel_error)
    b_clean = (
        sensor_coupling[:, 0 : len(mesh.vertices)] @ s.vert[0 : len(mesh.vertices)]
    )
    clean_error = np.linalg.norm(b_clean - field) / np.linalg.norm(field)
    print("Interference rejection relative error:", clean_error * 100, "%")
    if PLOT:
        plotter = fig = pv.Plotter()
        surf1 = plotter.add_mesh(
            triangle_mesh(*mesh.vertices.T, mesh.faces),
            scalars=s.vert[0 : len(mesh.vertices)],
            show_scalar_bar=False,
        )
        mesh_actor(surf1).mapper.interpolate_before_map = True
        mesh_actor(surf1).mapper.lookup_table.n_values = 16
        surf2 = plotter.add_mesh(
            triangle_mesh(*mesh2.vertices.T, mesh2.faces),
            scalars=s.vert[len(mesh.vertices) :],
            show_scalar_bar=False,
        )
        mesh_actor(surf2).mapper.dataset = mesh_actor(surf2).mapper.dataset.contour(10)
        currentloop.plot_loops(figure=fig)


#%%
# Interpolate MEG data to the sensor surface

from bfieldtools.utils import load_example_mesh

helmet = load_example_mesh("meg_helmet", process=False)
# Bring the surface roughly to the correct place
helmet.vertices[:, 2] -= 0.05

# Reset coupling by hand
c.B_coupling.reset()
plotter = pv.Plotter()
B_surf = np.sum(
    c.B_coupling(helmet.vertices) * helmet.vertex_normals[:, :, None], axis=1
)

if PLOT:
    fig = plotter
    sensors = plotter.add_mesh(
        vector_glyphs(*p.T, *n.T, mode="arrow"), show_scalar_bar=False
    )
    scalars = B_surf @ s
    surf = plotter.add_mesh(
        triangle_mesh(*helmet.vertices.T, helmet.faces),
        scalars=scalars,
        cmap="seismic",
        show_scalar_bar=False,
    )
    mesh_actor(surf).mapper.interpolate_before_map = True
    mesh_actor(surf).mapper.lookup_table.n_values = 15
    if SAVE_FIGURES:
        fig.screenshot(SAVE_DIR + "SUH_sensors_streamfunction.png", scale=4)


#%%
# Calculate magnetic field in volumetric grid

Nvol = 30
x = np.linspace(-0.125, 0.125, Nvol)
vol_points = np.array(np.meshgrid(x, x, x, indexing="ij")).reshape(3, -1).T

c.B_coupling.reset()
Bvol_coupling = c.B_coupling(vol_points, Nchunks=100, analytic=True)
s = StreamFunction(a, c)
# s = StreamFunction(a, c)
Bvol = Bvol_coupling @ s

#%%
# Now, plot the computed magnetic field using streamlines

if PLOT:
    from bfieldtools.mesh_calculus import gradient

    plotter = fig = pv.Plotter()
    plotter = surf = s.plot(False, figure=fig)
    mesh_actor(surf).mapper.interpolate_before_map = True
    mesh_actor(surf).mapper.lookup_table.n_values = 16
    vecs = structured_field(
        *vol_points.T.reshape(3, Nvol, Nvol, Nvol), *Bvol.T.reshape(3, Nvol, Nvol, Nvol)
    )
    vecnorm = vecs
    seed_points = mesh.vertices[mesh.faces].mean(axis=1) - 0.01 * mesh.face_normals
    seed_vals = c.basis @ c.inductance @ s
    seed_vals_grad = np.linalg.norm(gradient(seed_vals, c.mesh), axis=0)
    seed_vals = abs(seed_vals[mesh.faces].mean(axis=1)) ** 2
    seed_vals[seed_vals_grad > seed_vals_grad.max() / 1.8] = 0
    Npoints = 500
    seed_inds = np.random.choice(
        np.arange(len(seed_vals)), Npoints, False, seed_vals / seed_vals.sum()
    )
    seed_points = seed_points[seed_inds]
    streams = []
    for pi in seed_points:
        streamline = plotter.add_mesh(
            vecnorm.streamlines_from_source(
                pv.PolyData(pi.reshape(1, 3)),
                vectors="vectors",
                integration_direction="both",
                terminal_speed=1e-13,
                max_length=0.1,
            ),
            cmap="BuGn",
            show_scalar_bar=False,
        )
        mesh_actor(streamline).prop.render_lines_as_tubes = True
        mesh_actor(streamline).prop.line_width = 4.0
        streams.append(streamline)
    for streamline in streams:
        mesh_actor(streamline).mapper.lookup_table.apply_opacity(np.linspace(0, 1, 256))
        mesh_actor(streamline).mapper.scalar_range = (1e-13, 1e-12)
    sensors = plotter.add_mesh(
        vector_glyphs(*p.T, *n.T, mode="cylinder"), show_scalar_bar=False
    )
    mesh_actor(sensors).prop.color = (0.5, 0.5, 0.5)
    mesh_actor(sensors).mapper.scalar_visibility = False
    fig.camera.position = [0.637392177469018, 0.07644693029292644, -0.07183513804689762]
    fig.camera.focal_point = [
        -6.413459777832031e-05,
        0.01716560870409012,
        -0.0229007127850005,
    ]
    fig.camera.view_angle = 30.0
    fig.camera.up = [0.04390624852005244, 0.3114421192517664, 0.9492502555685007]
    fig.camera.clipping_range = [0.3366362817578398, 1.0281065506557443]
    camera_light = fig.renderer.lights[0]
    camera_light.intensity = 0.7
    if SAVE_FIGURES:
        fig.screenshot(SAVE_DIR + "SUH_streamlines_lateral.png", scale=4)
    fig.camera.position = [
        -6.413459777832031e-05,
        0.01716560870409012,
        0.6191735842078244,
    ]
    fig.camera.focal_point = [
        -6.413459777832031e-05,
        0.01716560870409012,
        -0.0229007127850005,
    ]
    fig.camera.view_angle = 30.0
    fig.camera.up = [0.0, 1.0, 0.0]
    fig.camera.clipping_range = [0.3381552363433513, 1.0261944997830243]
    if SAVE_FIGURES:
        fig.screenshot(SAVE_DIR + "SUH_streamlines_coronal.png", scale=4)


#%%
# Finally, plot the surface current density itself

plotter = fig = pv.Plotter()
grad_s = gradient(c.basis @ s, mesh, rotated=True)
q = plotter.add_mesh(
    vector_glyphs(*mesh.vertices[mesh.faces].mean(axis=1).T, *grad_s, mode="arrow"),
    cmap="viridis",
    show_scalar_bar=False,
)

plotter.add_mesh(
    triangle_mesh(*head.vertices.T, head.faces),
    color=(0.8, 0.8, 0.8),
    opacity=1.0,
    show_scalar_bar=False,
)

fig.camera.position = [0.4987072212753703, 0.06469079487766746, -0.0014732384935239248]
fig.camera.focal_point = [
    0.0018187984824180603,
    0.012344694641686624,
    -0.04367139294433087,
]
fig.camera.view_angle = 30.0
fig.camera.up = [-0.10720122151366927, 0.23975383168819672, 0.9648968848000314]
fig.camera.clipping_range = [0.28329092545021717, 0.7772019991936254]

if SAVE_FIGURES:
    fig.screenshot(SAVE_DIR + "SUH_surface_currents_lateral.png", scale=4)

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
