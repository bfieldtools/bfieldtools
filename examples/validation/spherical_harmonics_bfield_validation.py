"""
Spherical harmonics B-field computation validation
==================================================



"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import triangle_mesh, vector_glyphs

pv.set_plot_theme("document")


import numpy as np
import trimesh

from bfieldtools.mesh_magnetics import magnetic_field_coupling
from bfieldtools.mesh_conductor import MeshConductor

from bfieldtools.sphtools import compute_sphcoeffs_mesh
from bfieldtools import sphtools


from bfieldtools import resources

# Load simple plane mesh that is centered on the origin
file_obj = resources.files("bfieldtools", "example_meshes/10x10_plane.obj")
coilmesh = trimesh.load(file_obj, process=False)
coil = MeshConductor(mesh_obj=coilmesh)

coil.mesh.vertices += np.array([0, -1, 0])
weights = np.zeros(coilmesh.vertices.shape[0])
weights[coil.inner_vertices] = 1

test_points = coilmesh.vertices.copy()
test_points[:, 1] = 0

lmax = 9


sph_C = compute_sphcoeffs_mesh(coil.mesh, lmax)

alms = sph_C[0] @ weights
blms = sph_C[1] @ weights

alms = np.zeros_like(alms)

B0 = (magnetic_field_coupling(coilmesh, test_points) @ weights).T
B1 = sphtools.field(test_points, alms, blms, lmax).T

plotter = pv.Plotter()

s = plotter.add_mesh(
    triangle_mesh(*coilmesh.vertices.T, coilmesh.faces),
    scalars=weights,
    cmap="viridis",
    show_scalar_bar=False,
)
mesh_actor(s).mapper.dataset = mesh_actor(s).mapper.dataset.contour(10)
mesh_actor(s).prop.render_lines_as_tubes = True
mesh_actor(s).prop.line_width = 3.0

plotter.add_mesh(
    vector_glyphs(*test_points.T, *B0, scale_factor=5000000.0),
    color=(1, 0, 0),
    clim=[0, 2e-07],
    show_scalar_bar=False,
)
plotter.add_mesh(
    vector_glyphs(*test_points.T, *B1, scale_factor=5000000.0),
    color=(0, 0, 1),
    clim=[0, 2e-07],
    show_scalar_bar=False,
)
plotter.view_isometric()


print(
    "Relative RMS error", np.sqrt(np.mean((B1 - B0) ** 2)) / np.sqrt(np.mean((B0) ** 2))
)

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
