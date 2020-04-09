"""
Logo generation
=========================

This script generates the bfieldtools logo
"""

#
# import numpy as np
# import matplotlib.pyplot as plt
# import sys
# path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
# if path not in sys.path:
#    sys.path.insert(0,path)
#
# from bfieldtools.integrals import triangle_potential_dipole_linear
# from bfieldtools.integrals import omega
# from bfieldtools.utils import tri_normals_and_areas
# from bfieldtools.mesh_calculus import gradient
# from bfieldtools.mesh_magnetics import scalar_potential_coupling as compute_U
# from bfieldtools.mesh_magnetics import vector_potential_coupling as compute_A
# from bfieldtools.mesh_magnetics import magnetic_field_coupling as compute_C
from bfieldtools.mesh_magnetics import (
    magnetic_field_coupling_analytic as compute_C_analytic,
)

import trimesh
from mayavi import mlab
import numpy as np

#########################################################
# Test potential shape slightly above the surface

scalars = np.zeros(7)
scalars[0] = 1


def plot_element():
    x = np.sin(np.pi / 6)
    y = np.cos(np.pi / 6)
    points = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [x, y, 0],
            [-x, y, 0],
            [-1, 0, 0],
            [-x, -y, 0],
            [x, -y, 0],
        ]
    )

    tris = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 6], [0, 6, 1]])
    mesh = trimesh.Trimesh(points, tris)
    # Stream function
    mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap="Blues")
    # Stream lines
    s2 = mlab.triangular_mesh(*points.T, tris, scalars=scalars, colormap="Blues")
    s2.enable_contours = True
    s2.actor.mapper.scalar_range = np.array([0.0, 1.0])
    s2.actor.mapper.scalar_visibility = False
    s2.actor.property.render_lines_as_tubes = True
    s2.actor.property.line_width = 3.0
    s2.contour.number_of_contours = 6
    s2.actor.property.color = (1.0, 0.0, 0.0)

    return s2, mesh


s, mesh = plot_element()

#%%
# mlab.figure(bgcolor=(1,1,1))
#
# s, mesh = plot_element()
# points = np.array([[0.01, 1, 1],
#                   [0.01, 1, -1],
#                   [0.01, -1, -1],
#                   [0.01, -1, 1]])*2
# tris=np.array([[0,1,2], [2,3,0]])
# mesh2 = trimesh.Trimesh(points, tris)
# for ii in range(7):
#    mesh2 =mesh2.subdivide()
#
# U = compute_U(mesh, mesh2.vertices) @ scalars
#
# s3= mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=U, colormap='bwr')
# s3.enable_contours = True
# s3.contour.number_of_contours = 31
# s3.contour.filled_contours = True
# s3.contour.minimum_contour = -5.38e-07
# s3.contour.maximum_contour = 5.38e-07
##s3.actor.property.render_lines_as_tubes = True
# s3.actor.property.line_width = 3.0

#%%
# mlab.figure(bgcolor=(1,1,1))

# s, mesh = plot_element()
# points = np.array([[1, 1, -0.01],
#                   [1, -1, -0.01],
#                   [-1, -1, -0.01],
#                   [-1, 1, -0.01]])*2
# points[:,2] += 0.2
# tris=np.array([[0,1,2], [2,3,0]])
# mesh3 = trimesh.Trimesh(points, tris)
# for ii in range(3):
#    mesh3 =mesh3.subdivide()
# A = compute_A(mesh, mesh3.vertices) @ scalars
# vectors = mlab.quiver3d(*mesh3.vertices.T, *A, mode='arrow', color=(0,0,1))
# vectors.glyph.glyph_source.glyph_position = 'center'
#%%
# mlab.figure(bgcolor=(1,1,1))
# s, mesh = plot_element()
# points = np.array([[0.001, 1, 1],
#                   [0.001, 1, -1],
#                   [0.001, -1, -1],
#                   [0.001, -1, 1]])*2 + 0.001
p1 = np.array([-0.1, -1, 0]) * 2 + 0.0001
p2 = np.array([0.1, 1, 1]) * 2 + 0.0001
X, Y, Z = np.meshgrid(*np.linspace(p1, p2, 50).T, indexing="ij")
R = np.array([X.flatten(), Y.flatten(), Z.flatten()])
# B0 = compute_C(mesh, R.T) @ scalars
B1 = compute_C_analytic(mesh, R.T) @ scalars
B1[0] = 0
vecfield = mlab.pipeline.vector_field(X, Y, Z, *B1.T.reshape(3, 50, 50, 50))
streamline = mlab.pipeline.streamline(vecfield)
streamline.seed.widget = streamline.seed.widget_list[1]
# streamline.seed.widget = <tvtk.tvtk_classes.line_widget.LineWidget object at 0x7f96be04f678>
streamline.stream_tracer.start_position = np.array([0.0, 0.0, 0.0])
streamline.stream_tracer.integration_direction = "both"
streamline.seed.widget.resolution = 15
streamline.seed.widget.point1 = np.array([0.01, -0.6, 0.055])
streamline.seed.widget.point2 = np.array([0.01, 0.6, 0.05523521])
streamline.actor.property.render_lines_as_tubes = True
streamline.actor.property.line_width = 3.0
streamline.actor.property.color = (0, 0, 0)
streamline.update_streamlines = 0

# The other half
p1 = np.array([-0.1, -1, -1]) * 2 - 0.0001
p2 = np.array([0.1, 1, 0]) * 2 - 0.0001
X, Y, Z = np.meshgrid(*np.linspace(p1, p2, 50).T, indexing="ij")
R = np.array([X.flatten(), Y.flatten(), Z.flatten()])
# B0 = compute_C(mesh, R.T) @ scalars
B1 = compute_C_analytic(mesh, R.T) @ scalars
B1[0] = 0
vecfield = mlab.pipeline.vector_field(X, Y, Z, *B1.T.reshape(3, 50, 50, 50))
streamline = mlab.pipeline.streamline(vecfield)
streamline.seed.widget = streamline.seed.widget_list[1]
# streamline.seed.widget = <tvtk.tvtk_classes.line_widget.LineWidget object at 0x7f96be04f678>
streamline.stream_tracer.start_position = np.array([0.0, 0.0, 0.0])
streamline.stream_tracer.integration_direction = "both"
streamline.seed.widget.resolution = 15
streamline.seed.widget.point1 = np.array([0.01, -0.6, -0.055])
streamline.seed.widget.point2 = np.array([0.01, 0.6, -0.05523521])
streamline.actor.property.render_lines_as_tubes = True
streamline.actor.property.line_width = 3.0
streamline.actor.property.color = (0, 0, 0)
streamline.update_streamlines = 0


# vectors = mlab.quiver3d(*R, *B1.T, mode='arrow', color=(1,0,1))
# vectors.glyph.glyph_source.glyph_position = 'center'
# vectors.actor.property.render_lines_as_tubes = True
# vectors.actor.property.line_width = 3.0
