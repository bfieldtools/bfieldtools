"""
Figure 2: Surface harmonics modes
==================================

First 15 surface harmonics of a triangle mesh. The index (and spatial frequency) increases from left to right, row by row. The hole and outer boundary have zero tangential derivative. The mesh discretization is shown in the magnified inset to the right.

"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor

pv.set_plot_theme("document")

from bfieldtools.mesh_conductor import MeshConductor
from bfieldtools import resources

import numpy as np

SAVE = False
SAVE_DIR = "./Surface harmonics/"


c = MeshConductor(
    mesh_file=resources.files("bfieldtools", "example_meshes/curved_surf_w_hole.stl"),
    process=True,
    basis_name="suh",
    N_suh=15,
    fix_normals=True,
)


T_x = 1.5 * np.pi / 2
T_z = -1.02 * np.pi
rotmat = np.array(
    [
        [np.cos(T_z), -np.sin(T_z), 0, 0],
        [np.sin(T_z), np.cos(T_z), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
) @ np.array(
    [
        [1, 0, 0, 0],
        [0, np.cos(T_x), -np.sin(T_x), 0],
        [0, np.sin(T_x), np.cos(T_x), 0],
        [0, 0, 0, 1],
    ]
)


c.mesh.apply_transform(rotmat)


plotter = f = pv.Plotter(window_size=(750, 900))

surfaces = c.suh_basis.plot(
    Nfuncs=c.basis.shape[1], Ncols=2, dist=0.1, colormap="RdBu", figure=f, ncolors=256
)

surfaces = [a for a in f.actors.values() if isinstance(a, pv.Actor)]
f.view_xy()
1 + 1
#
for (idx, surface) in enumerate(surfaces):
    mesh_actor(surface).mapper.interpolate_before_map = True
    mesh_actor(surface).mapper.scalar_range = np.array(
        [
            -np.max(np.abs(c.suh_basis.basis[:, idx])),
            np.max(np.abs(c.suh_basis.basis[:, idx])),
        ]
    )


f.camera.zoom(1.2)

if SAVE:
    f.screenshot(SAVE_DIR + "surface_harmonics.png", scale=4)
    plotter.close()

plotter = f = pv.Plotter(window_size=(750, 600))

plotter = c.plot_mesh(style="wireframe", figure=f)
plotter = c.plot_mesh(opacity=0.2, figure=f)
#
f.camera.zoom(1.5)
f.view_xy()

if SAVE:
    f.screenshot(SAVE_DIR + "suhmesh.png", scale=6)
    plotter.close()

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
