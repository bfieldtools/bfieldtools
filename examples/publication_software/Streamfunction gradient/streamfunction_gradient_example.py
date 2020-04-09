#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 12:44:50 2020

@author: Rasmus Zetter
"""

from bfieldtools.conductor import Conductor, StreamFunction
import pkg_resources
from bfieldtools.mesh_calculus import gradient
import numpy as np

c = Conductor(
    mesh_file=pkg_resources.resource_filename(
        "bfieldtools", "example_meshes/curved_surf_w_hole.stl"
    ),
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


from mayavi import mlab
from mayavi.api import Engine

engine = Engine()
engine.start()

f = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 700))
s.plot(figure=f, ncolors=256)
c.plot_mesh(representation="wireframe", figure=f)


j = gradient(s.vert, c.mesh, rotated=True)

Len = np.log(np.linalg.norm(j, axis=0))

mlab.quiver3d(
    *c.mesh.triangles_center.T, *j, mode="arrow", colormap="Greens", scalars=Len
)

vectors = engine.scenes[0].children[2].children[0].children[0]
vectors.glyph.glyph.scale_mode = "scale_by_scalar"
vectors.glyph.glyph.scale_factor = 0.6
f.scene.z_plus_view()


module_manager2 = engine.scenes[0].children[2].children[0]
module_manager2.scalar_lut_manager.scalar_bar_representation.maximum_size = np.array(
    [100000, 100000]
)
module_manager2.scalar_lut_manager.scalar_bar_representation.minimum_size = np.array(
    [1, 1]
)
module_manager2.scalar_lut_manager.scalar_bar_representation.position = np.array(
    [0.82, 0.1]
)
module_manager2.scalar_lut_manager.scalar_bar_representation.position2 = np.array(
    [0.17, 0.8]
)
module_manager2.scalar_lut_manager.show_scalar_bar = True
module_manager2.scalar_lut_manager.show_legend = True
module_manager2.scalar_lut_manager.scalar_bar.height = 0.8
module_manager2.scalar_lut_manager.scalar_bar.position = np.array([0.82, 0.1])
module_manager2.scalar_lut_manager.scalar_bar.position2 = np.array([0.17, 0.8])
module_manager2.scalar_lut_manager.scalar_bar.width = 0.17
scene = engine.scenes[0]
scene.scene.camera.position = [
    -0.3696892487983681,
    0.2840788710848503,
    3.701830880912346,
]
scene.scene.camera.focal_point = [
    -0.3696892487983681,
    0.2840788710848503,
    0.8575533408480627,
]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [2.5164461179149695, 3.263810326333801]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
module_manager2.scalar_lut_manager.scalar_bar.position = np.array([0.82, 0.1])
module_manager2.scalar_lut_manager.scalar_bar.position2 = np.array([0.17, 0.8])
module_manager2.scalar_lut_manager.scalar_bar.number_of_labels = 0
module_manager2.scalar_lut_manager.scalar_bar.position = np.array([0.82, 0.1])
module_manager2.scalar_lut_manager.scalar_bar.position2 = np.array([0.17, 0.8])
module_manager2.scalar_lut_manager.number_of_labels = 0
module_manager2.scalar_lut_manager.scalar_bar.position = np.array([0.82, 0.1])
module_manager2.scalar_lut_manager.scalar_bar.position2 = np.array([0.17, 0.8])
module_manager2.scalar_lut_manager.scalar_bar.maximum_number_of_colors = 8
module_manager2.scalar_lut_manager.scalar_bar.position = np.array([0.82, 0.1])
module_manager2.scalar_lut_manager.scalar_bar.position2 = np.array([0.17, 0.8])
module_manager2.scalar_lut_manager.number_of_colors = 8
module_manager2.scalar_lut_manager.use_default_name = False
module_manager2.scalar_lut_manager.scalar_bar.position = np.array([0.82, 0.1])
module_manager2.scalar_lut_manager.scalar_bar.position2 = np.array([0.17, 0.8])
module_manager2.scalar_lut_manager.scalar_bar.title = "Current density"
module_manager2.scalar_lut_manager.scalar_bar.position = np.array([0.82, 0.1])
module_manager2.scalar_lut_manager.scalar_bar.position2 = np.array([0.17, 0.8])
module_manager2.scalar_lut_manager.data_name = "Current density"
module_manager2.scalar_lut_manager.label_text_property.shadow_offset = np.array([1, -1])
module_manager2.scalar_lut_manager.label_text_property.italic = False
module_manager2.scalar_lut_manager.label_text_property.shadow_offset = np.array([1, -1])
module_manager2.scalar_lut_manager.label_text_property.color = (0.0, 0.0, 0.0)
module_manager2.scalar_lut_manager.title_text_property.shadow_offset = np.array([1, -1])
module_manager2.scalar_lut_manager.title_text_property.italic = False
module_manager2.scalar_lut_manager.title_text_property.shadow_offset = np.array([1, -1])
module_manager2.scalar_lut_manager.title_text_property.color = (0.0, 0.0, 0.0)
module_manager2.scalar_lut_manager.title_text_property.shadow_offset = np.array([1, -1])
module_manager2.scalar_lut_manager.title_text_property.bold = False


mlab.savefig(
    "/l/bfieldtools/examples/publication_software/streamfunction_gradient.png",
    figure=f,
    magnification=4,
)
