"""
Integrals testing
==================================================

"""

import numpy as np

import matplotlib.pyplot as plt

import trimesh
from mayavi import mlab

#%%
# Test potential shape slightly above the surface
#%%
x = np.sin(np.pi / 6)
y = np.cos(np.pi / 6)
points = (
    np.array(
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
    * 2
)

tris = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 6], [0, 6, 1]])
mesh = trimesh.Trimesh(points, tris)
scalars = np.zeros(7)
scalars[0] = 1

#%%
# Linear dipole density
#%%

# Sign ok
points = np.array([[0.1, 1, 1], [0.1, 1, -1], [0.1, -1, -1], [0.1, -1, 1]]) * 2
# points = np.roll(points, 2, 1)
tris = np.array([[0, 1, 2], [2, 3, 0]])
mesh2 = trimesh.Trimesh(points, tris)
for ii in range(7):
    mesh2 = mesh2.subdivide()

from bfieldtools.legacy.integrals import triangle_potential_dipole_linear as t1
from bfieldtools.integrals import triangle_potential_dipole_linear as t2

RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
p1 = t1(RR, mesh.face_normals, mesh.area_faces, planar=False)
p2 = t2(RR, mesh.face_normals, mesh.area_faces)

assert np.allclose(p1, p2)


mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p1[:, :, 0].sum(axis=1))
mlab.colorbar()
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p2[:, :, 0].sum(axis=1))
mlab.colorbar()
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5), size=(800, 800))
mlab.triangular_mesh(
    *mesh2.vertices.T, mesh2.faces, scalars=(p1 - p2)[:, :, 0].sum(axis=1)
)
mlab.colorbar()

#%%
#
points = np.zeros((100, 3))
points[:, 2] = np.linspace(-1, 1, 100)
from bfieldtools.legacy.integrals import omega as omega1
from bfieldtools.integrals import omega as omega2

RR = points[:, None, None, :] - mesh.vertices[None, mesh.faces]
o1 = omega1(RR).sum(axis=1)
o2 = omega2(RR).sum(axis=1)

assert np.allclose(o1, -o2)

plt.plot(o1)
plt.plot(o2)
mlab.plot3d(*points.T, points[:, 2], colormap="seismic")
mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

#%%
# Plot x_i

from bfieldtools.integrals import x_distance

RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
xdist = x_distance(RR, mesh.face_normals)
mlab.triangular_mesh(
    *mesh2.vertices.T,
    mesh2.faces,
    scalars=xdist[:, 1, 0],
    vmin=-1,
    vmax=1,
    colormap="seismic"
)
mlab.colorbar()
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

#%%
# Uniform charge density
from bfieldtools.legacy.integrals import triangle_potential_uniform as u1
from bfieldtools.integrals import triangle_potential_uniform as u2

RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
p1 = u1(RR, mesh.face_normals, planar=False)
p2 = u2(RR, mesh.face_normals, planar=False)

assert np.allclose(p1, p2)


mlab.figure("uniform charge density (old)")
mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p1.sum(axis=1))
mlab.colorbar()
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)
mlab.figure("uniform charge density (new)")
mlab.triangular_mesh(*mesh2.vertices.T, mesh2.faces, scalars=p2.sum(axis=1))
mlab.colorbar()
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

#%%
#
from bfieldtools.integrals import d_distance

RR = mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
ddist = d_distance(RR, mesh.face_normals)
mlab.figure("d distance")
mlab.triangular_mesh(
    *mesh2.vertices.T,
    mesh2.faces,
    scalars=ddist[:, 0],
    vmin=-1,
    vmax=1,
    colormap="seismic"
)
mlab.colorbar()
mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, representation="wireframe")
mlab.quiver3d(*mesh.triangles_center.T, *mesh.face_normals.T)

#%%
#
from bfieldtools.legacy.mesh_magnetics import (
    magnetic_field_coupling_analytic as magnetic_field_coupling_analytic_old,
)
from bfieldtools.mesh_magnetics import magnetic_field_coupling_analytic

b1 = magnetic_field_coupling_analytic_old(mesh, mesh2.vertices)
b2 = magnetic_field_coupling_analytic(mesh, mesh2.vertices)

assert np.allclose(b1, b2)

mlab.figure("b field")
mlab.quiver3d(*mesh2.vertices.T, *b1[:, :, 0].T)
mlab.quiver3d(*mesh2.vertices.T, *b2[:, :, 0].T)

#%%
# Gammma
from bfieldtools.legacy.integrals import gamma0 as g1
from bfieldtools.integrals import gamma0 as g2

# RR =  mesh2.vertices[:, None, None, :] - mesh.vertices[None, mesh.faces]
t = np.linspace(-1.5, 1.5)
points = (
    t[:, None] * mesh.vertices[mesh.faces][0][0]
    + (1 - t)[:, None] * mesh.vertices[mesh.faces][0][1]
)


R = points[:, None, None, :] - mesh.vertices[None, mesh.faces]
p1 = g1(R, symmetrize=True)
p2 = g2(R, symmetrize=True)

assert np.allclose(p1, p2)

plt.figure()
plt.plot(p1[:, 0, :])
plt.figure()
plt.plot(p2[:, 0, :])
