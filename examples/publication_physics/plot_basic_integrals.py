"""
Analytic integrals as scalar fields
========================================================

"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh, vector_glyphs, polyline

pv.set_plot_theme("document")
plotter = pv.Plotter()


import numpy as np
import matplotlib.pyplot as plt
import trimesh


from bfieldtools.integrals import triangle_potential_uniform
from bfieldtools.integrals import triangle_potential_dipole_linear
from bfieldtools.integrals import gamma0
from bfieldtools.integrals import omega
from bfieldtools.utils import tri_normals_and_areas

#%%
#%% Create there orthogonal planes
points = np.array([[0, 0, 0], [1, 0.5, 0], [0, 1, 0]]) + 0.15

points[:, 2] += 0.1

tris = np.array([[0, 1, 2]])
p_tris = points[tris]
mesh = trimesh.Trimesh(points, tris, process=False)

# Evaluation points
points2 = np.array([[-0.01, 1, 1], [-0.01, 1, 0], [-0.01, 0, 0], [-0.01, 0, 1]]) * 1.2
tris2 = np.array([[0, 1, 2], [2, 3, 0]])
mesh2 = trimesh.Trimesh(points2, tris2)
for ii in range(7):
    mesh2 = mesh2.subdivide()

points3 = np.array([[1, 1, -0.01], [1, 0, -0.01], [0, 0, -0.01], [0, 1, -0.01]]) * 1.2
mesh3 = trimesh.Trimesh(points3, tris2)
for ii in range(7):
    mesh3 = mesh3.subdivide()

points4 = np.array([[1, -0.01, 1], [1, -0.01, 0], [0, -0.01, 0], [0, -0.01, 1]]) * 1.2
mesh4 = trimesh.Trimesh(points4, tris2)
for ii in range(7):
    mesh4 = mesh4.subdivide()


def plot_axes():
    axes_points = np.zeros((3, 3))
    axes = np.eye(3) * 0.3
    vectors = plotter.add_mesh(
        vector_glyphs(*axes_points, *axes, mode="arrow"),
        color=(0, 0, 0),
        show_scalar_bar=False,
    )


# Difference vectors
RR2 = mesh2.vertices[:, None, None, :] - p_tris[None, :, :, :]
RR3 = mesh3.vertices[:, None, None, :] - p_tris[None, :, :, :]
RR4 = mesh4.vertices[:, None, None, :] - p_tris[None, :, :, :]

tn, ta = tri_normals_and_areas(points, tris)

#%% Plot potentials on the planes and the respective sources on the triangle
for (ii, func) in enumerate(
    (triangle_potential_uniform, triangle_potential_dipole_linear, gamma0, omega)
):
    plotter = pv.Plotter()
    print(func)
    if ii == 0:
        pot2 = func(RR2, tn)[:, 0]
        pot3 = func(RR3, tn)[:, 0]
        pot4 = func(RR4, tn)[:, 0]
        plotter.add_mesh(
            triangle_mesh(*mesh.vertices.T, mesh.faces),
            color=(0.5, 0.5, 0.5),
            opacity=0.7,
            show_scalar_bar=False,
        )
    if ii == 1:
        pot2 = func(RR2, tn, ta)[:, 0, 1]
        pot3 = func(RR3, tn, ta)[:, 0, 1]
        pot4 = func(RR4, tn, ta)[:, 0, 1]
        plotter.add_mesh(
            triangle_mesh(*mesh.vertices.T, mesh.faces),
            color=(0.5, 0.5, 0.5),
            opacity=0.2,
            show_scalar_bar=False,
        )
        meshq = mesh.copy()
        for jj in range(4):
            meshq = meshq.subdivide()
        u = np.zeros(meshq.vertices.shape)
        r = meshq.vertices
        r2 = np.zeros(r.shape + (3,))
        r2[:, 1] = r
        r2[:, 0] = mesh.vertices[0]
        r2[:, 2] = mesh.vertices[2]
        u[:, 2] = np.linalg.det(r2) / np.linalg.det(mesh.vertices)
        q = plotter.add_mesh(
            vector_glyphs(*r.T, *u.T, mode="arrow"), cmap="gray", show_scalar_bar=False
        )
    if ii == 2:
        pot2 = func(RR2)[:, 0, 2]
        pot3 = func(RR3)[:, 0, 2]
        pot4 = func(RR4)[:, 0, 2]
        plotter.add_mesh(
            polyline(*points[0:2].T).tube(radius=0.02),
            color=(0.5, 0.5, 0.5),
            opacity=0.5,
            show_scalar_bar=False,
        )
    if ii == 3:
        pot2 = func(RR2)[:, 0]
        pot3 = func(RR3)[:, 0]
        pot4 = func(RR4)[:, 0]
        plotter.add_mesh(
            triangle_mesh(*mesh.vertices.T, mesh.faces),
            color=(0.5, 0.5, 0.5),
            opacity=0.5,
            show_scalar_bar=False,
        )
        u = np.zeros(meshq.vertices.shape)
        u[:, 2] = 1
        r = meshq.vertices
        q = plotter.add_mesh(
            vector_glyphs(*r.T, *u.T, mode="arrow"), cmap="gray", show_scalar_bar=False
        )
    plot_axes()
    M = max(max(abs(pot2)), max(abs(pot3)), max(abs(pot4)))
    for (m, p) in zip((mesh2, mesh3, mesh4), (pot2, pot3, pot4)):
        s = plotter.add_mesh(
            triangle_mesh(*m.vertices.T, m.faces),
            scalars=p,
            cmap="seismic",
            clim=[-M, M],
            show_scalar_bar=False,
        )
        mesh_actor(s).mapper.interpolate_before_map = True
        mesh_actor(s).mapper.lookup_table.n_values = 32
    plotter.camera.Pitch(-2)
    plotter.camera.Azimuth(-25)
    plotter.camera.Elevation(-10)
    plotter.camera.OrthogonalizeViewUp()
    plotter.camera.Pitch(1)
    plotter.camera.OrthogonalizeViewUp()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
