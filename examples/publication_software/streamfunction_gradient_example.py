"""
Figure 1: Stream function example
==================================

An example stream function (red--blue colormap) and its rotated gradient, i.e. the surface current density (arrows; green colormap) on a surface mesh with a hole in it. The surface normal is oriented up towards the reader.
"""

import pyvista as pv
from bfieldtools._example_plots import mesh_actor
from bfieldtools._example_plots import vector_glyphs

pv.set_plot_theme("document")


from bfieldtools.mesh_conductor import MeshConductor, StreamFunction
from bfieldtools import resources
from bfieldtools.mesh_calculus import gradient
import numpy as np

SAVE = False
SAVE_DIR = "./Streamfunction gradient/"


c = MeshConductor(
    mesh_file=resources.files("bfieldtools", "example_meshes/curved_surf_w_hole.stl"),
    process=True,
    basis_name="suh",
    N_suh=10,
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

s = np.zeros((c.basis.shape[1],))
s[2] += 1
# s[63] += 2

s = StreamFunction(s, c)


plotter = f = pv.Plotter(window_size=(800, 700))
plotter = s.plot(figure=f, ncolors=256)
plotter = c.plot_mesh(style="wireframe", figure=f)


j = gradient(s.vert, c.mesh, rotated=True)

Len = np.log(np.linalg.norm(j, axis=0))

vectors = plotter.add_mesh(
    vector_glyphs(*c.mesh.triangles_center.T, *j, mode="arrow", scalars=Len),
    cmap="Greens",
    show_scalar_bar=False,
)

f.view_xy()

#
mesh_actor(vectors).mapper.lookup_table.n_values = 8
plotter.add_scalar_bar(title="Current density", n_labels=0)
f.camera.position = [-0.3696892487983681, 0.2840788710848503, 3.701830880912346]
f.camera.focal_point = [-0.3696892487983681, 0.2840788710848503, 0.8575533408480627]
f.camera.view_angle = 30.0
f.camera.up = [0.0, 1.0, 0.0]
f.camera.clipping_range = [2.5164461179149695, 3.263810326333801]
f.render()

if SAVE:
    f.screenshot(SAVE_DIR + "streamfunction_gradient.png", scale=4)

# Display the assembled scene (Sphinx-Gallery captures it during docs builds).
if "plotter" in globals() and not plotter._closed:
    plotter.show()
