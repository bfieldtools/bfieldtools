"""
Visualize SUH components on 3 surfaces
=====================================================
"""

import pyvista as pv
from bfieldtools._example_plots import triangle_mesh

pv.set_plot_theme("document")


import numpy as np
import matplotlib.pyplot as plt
import trimesh

from bfieldtools.suhtools import SuhBasis
from bfieldtools.utils import find_mesh_boundaries

from trimesh.creation import icosphere
from bfieldtools.utils import load_example_mesh

# Import meshes

# Sphere
sphere = icosphere(4, 0.1)

# Plane mesh centered on the origin
plane = load_example_mesh("10x10_plane_hires")
scaling_factor = 0.02
plane.apply_scale(scaling_factor)
# Rotate to x-plane
t = np.eye(4)
t[1:3, 1:3] = np.array([[0, 1], [-1, 0]])
plane.apply_transform(t)

# Bunny
bunny = load_example_mesh("bunny_repaired")
bunny.vertices -= bunny.vertices.mean(axis=0)

plotter = pv.Plotter()
for mesh in (sphere, plane, bunny):
    s = plotter.add_mesh(
        triangle_mesh(*mesh.vertices.T, mesh.faces),
        color=(0.5, 0.5, 0.5),
        opacity=0.2,
        show_scalar_bar=False,
    )
plotter.view_xy()

Nc = 20
basis_sphere = SuhBasis(sphere, Nc)
basis_plane = SuhBasis(plane, Nc)
basis_bunny = SuhBasis(bunny, Nc)


A = [mesh.area for mesh in (bunny, sphere, plane)]

#%%

Nfuncs = [0, 1, 2, 3, 4, 5]
kwargs = {"colormap": "RdBu", "ncolors": 15}
plt.figure(figsize=(3.5, 2.5))
for (i, b) in enumerate((basis_bunny, basis_sphere, basis_plane)):
    plotter = fig = pv.Plotter(window_size=(550, 150))
    s = b.plot(Nfuncs, 0.1, Ncols=6, figure=fig, **kwargs)
    plotter.parallel_projection = True
    plotter.view_xy()
    if i == 0:
        plotter.camera.parallel_scale = 0.1
    else:
        plotter.camera.parallel_scale = 0.13
    plt.plot(np.sqrt(b.eigenvals * A[i]), ".-")

plt.legend(("bunny", "sphere", "square"), loc="lower right")
plt.xlabel("component index $n$")
plt.ylabel("$\sqrt{A}k_n$")
plt.tight_layout()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
