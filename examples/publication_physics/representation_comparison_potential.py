"""

Compare SUH and SPH basis functions for the magnetic scalar potential
======================================================================

"""

import pyvista as pv

pv.set_plot_theme("document")


import numpy as np
import matplotlib.pyplot as plt
import trimesh
from bfieldtools import resources


from bfieldtools.mesh_magnetics import magnetic_field_coupling
from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic
from bfieldtools.mesh_magnetics import scalar_potential_coupling
from bfieldtools.sphtools import compute_sphcoeffs_mesh, basis_fields
from bfieldtools.suhtools import SuhBasis
from bfieldtools.utils import load_example_mesh


mesh = load_example_mesh("bunny_repaired", process=True)

mesh.vertices -= mesh.vertices.mean(axis=0)

mesh_field = mesh.copy()
mesh_field.vertices += 0.001 * mesh_field.vertex_normals
# mesh_field = trimesh.smoothing.filter_laplacian(mesh_field, iterations=1)

# Ca, Cb = basis_fields(mesh_field.vertices, 4)

# bsuh = SuhBasis(mesh, 25)
# Csuh = magnetic_field_coupling_analytic(mesh, mesh_field.vertices) @ bsuh.basis

#%%
# def plot_basis_fields(C, comps):
#     d = 0.17
#     i = 0
#     j = 0
#     for n in comps:
#         p = 1.05 * mesh_field.vertices.copy()
#         p2 = mesh_field.vertices.copy()
#         #        p[:,1] -= i*d
#         #        p2[:,1] -= i*d
#         p[:, 0] += i * d
#         p2[:, 0] += i * d
#         m = np.max(np.linalg.norm(C[:, :, n], axis=0))
#             *p.T, *C[:, :, n].T, mode="arrow", colormap="Greys", vmin=0, vmax=m
#         )
#         #        m = np.max(abs((C[:,:,n].T*mesh_field.vertex_normals.T).sum(axis=0)))
#         #                             scalars=(C[:,:,n].T*mesh_field.vertex_normals.T).sum(axis=0),
#         #                             colormap='seismic', vmin=-m, vmax=m, opacity=0.7)
#         m = np.max(abs((C[:, :, n].T * mesh_field.vertex_normals.T).sum(axis=0)))
#             *p2.T,
#             mesh.faces,
#             scalars=(C[:, :, n].T * mesh_field.vertex_normals.T).sum(axis=0),
#             colormap="bwr",
#             vmin=-m,
#             vmax=m
#         )
#         i += 1


# comps = [0, 4, 10, 15]
# plot_basis_fields(Ca, comps)


# plot_basis_fields(Csuh, comps)


#%% test representation
from bfieldtools.mesh_magnetics import scalar_potential_coupling
from bfieldtools.sphtools import basis_potentials

scaling_factor = 0.02
# Load simple plane mesh that is centered on the origin
plane = load_example_mesh("10x10_plane_hires")
plane.apply_scale(scaling_factor)
# Rotate to x-plane
t = np.eye(4)
t[1:3, 1:3] = np.array([[0, 1], [-1, 0]])
plane.apply_transform(t)
# Subdivide face close to the mesh
from trimesh.proximity import signed_distance

dd = signed_distance(mesh, plane.triangles_center)
plane = plane.subdivide(np.flatnonzero(abs(dd) < 0.005))
dd = signed_distance(mesh, plane.triangles_center)
plane = plane.subdivide(np.flatnonzero(abs(dd) < 0.002))

UB = scalar_potential_coupling(mesh, plane.vertices, multiply_coeff=False)
mask = np.sum(UB, axis=1)
mask = mask > -0.5

Ua, Ub = basis_potentials(plane.vertices, 6)
bsuh = SuhBasis(mesh, 48)
# UB = scalar_potential_coupling(mesh, plane.vertices)
Usuh = UB @ bsuh.basis

sphere = trimesh.creation.icosphere(radius=0.02)
Ua_mesh, Ub_mesh = basis_potentials(sphere.vertices, 5)
UB_mesh = scalar_potential_coupling(mesh, mesh_field.vertices)
Usuh_mesh = UB_mesh @ bsuh.basis


#%%
comps = [0, 5, 13, 17]
d = 0.22
from bfieldtools.viz import plot_data_on_vertices
from bfieldtools.viz import plot_mesh

# Plot suh
i = 0
plotter = fig = pv.Plotter(window_size=(600, 190))
for n in comps:
    p = plane.copy()
    p.vertices[:, 0] += i * d
    p2 = mesh_field.copy()
    p2.vertices[:, 0] += i * d
    scalars = Usuh[:, n].copy()
    scalars[~mask] = 0
    scalars2 = bsuh.basis[:, n]
    m = np.max(abs(scalars))
    m2 = np.max(abs(scalars2))
    plotter = plot_data_on_vertices(
        p, scalars, figure=fig, ncolors=15, vmax=m, colormap="bwr"
    )
    plotter = plot_data_on_vertices(
        p2, scalars2, figure=fig, ncolors=15, vmax=m2, colormap="BrBG"
    )
    i += 1
fig.parallel_projection = True
fig.view_xy()
fig.camera.parallel_scale = 0.11

#%% Plot SPH
i = 0
plotter = fig = pv.Plotter(window_size=(600, 190))
for n in comps:
    p = plane.copy()
    p.vertices[:, 0] += i * d
    p2 = sphere.copy()
    p2.vertices[:, 0] += i * d
    scalars = Ua[:, n].copy()
    scalars[~mask] = 0
    scalars2 = Ua_mesh[:, n]
    m = np.max(abs(scalars))
    m2 = np.max(abs(scalars2))
    plotter = plot_data_on_vertices(
        p, scalars, figure=fig, ncolors=15, vmax=m, colormap="bwr"
    )
    plotter = plot_data_on_vertices(
        p2, scalars2, figure=fig, ncolors=15, vmax=m2, colormap="BrBG"
    )
    p3 = mesh.copy()
    p3.vertices[:, 0] += i * d
    plotter = plot_mesh(p3, opacity=0.3, figure=fig)
    i += 1
fig.parallel_projection = True
fig.view_xy()
fig.camera.parallel_scale = 0.11

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
