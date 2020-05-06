# -*- coding: utf-8 -*-
"""
Flatten mesh using conformal mapping
===========================================

Map 3D mesh to a 2D (complex) plane with angle-preserving (conformal) mapping

Based on these course notes
https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
section 7.4.

"""

import numpy as np
from bfieldtools.mesh_calculus import laplacian_matrix, mass_matrix
from bfieldtools.utils import find_mesh_boundaries, load_example_mesh
from scipy.sparse.linalg import eigsh
from scipy.sparse import csc_matrix
import trimesh


def A_matrix_complex(mesh):
    """
    Area matrix for complex parametrization
    See this course:
    https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
    section 7.4.

    Parameters
    ----------
    mesh : trimesh

    Returns
    -------
    A : area matrix
        u.T.conj() @ A @ u is the area of the mesh
        where u is the complex coordinates of the vertices 
        (complex parameterization)

    """
    p = mesh.outline().entities[0]
    b = p.points
    N = mesh.vertices.shape[0]
    ii = np.hstack([b, np.roll(b, 1)])
    jj = np.hstack([np.roll(b, 1), b])
    vv = np.ones(len(b))
    data = -0.25 * 1j * np.hstack([vv, -vv])
    A = csc_matrix((data, (ii, jj)), shape=(N, N))

    return A


def flatten_mesh(mesh, _lambda=1.0):
    """ Flatten the mesh, return uv coordinates and the mesh in 2D  

    Parameters
    ----------
    mesh : Trimesh object 
        must have boundary
    _lambda : int <= 1.0
        parameter for trading of area-distortion/angle-preservation. 
        The default is 1.0

    Returns
    -------
    u : array
        first coordinate of the paramterization
    v : array
        second coordinate of the paramterization
    mesh2d : Trimesh object with coordinates (u,v,0)

    lambda <= 1.0
    if lambda == 1.0 => conformal mapping
    if lambda == 0.5 =>  not conformal but less area distortion
    if lambda > 1 (e.g. 1.01-1.1) weird folding effects
    """

    N = mesh.vertices.shape[0]

    L = laplacian_matrix(mesh)
    M = mass_matrix(mesh)
    Ac = A_matrix_complex(mesh)
    vals, uv = eigsh(-0.5 * L.T - _lambda * Ac, 6, M, which="LM", sigma=0)

    # Coordinates with initial phase
    u = uv[:, 1].real
    v = uv[:, 1].imag

    # Determine "phase" by matching the uv coordinate function with mesh coordinates
    theta = np.linspace(0, 2 * np.pi, 50)
    yy = np.imag(np.exp(1j * theta)[:, None] * uv[:, 1])
    # plt.plot(np.sum(mesh.vertices[:,0]*xx, axis=1))
    ii = np.argmax(np.sum(mesh.vertices[:, 1] * yy, axis=1))

    theta = theta[ii]
    u = np.real(np.exp(1j * theta) * uv[:, 1])
    v = np.imag(np.exp(1j * theta) * uv[:, 1])

    mesh2d = trimesh.Trimesh(np.array([u, v, 0 * u]).T, mesh.faces, process=False)
    return u, v, mesh2d


# Map points from 2D to 3D or vice versa


def mesh2plane(points3d, mesh, u, v):
    c, d, f = trimesh.proximity.closest_point(mesh, points3d)
    tris = mesh.vertices[mesh.faces[f]]
    barys = trimesh.triangles.points_to_barycentric(tris, c)
    print(barys)
    p1 = np.sum(u[mesh.faces[f]] * barys, axis=1)
    p2 = np.sum(v[mesh.faces[f]] * barys, axis=1)
    return np.array([p1, p2]).T


def plane2mesh(points2d, mesh, u, v):
    mesh2d = trimesh.Trimesh(np.array([u, v, 0 * u]).T, mesh.faces)
    c, d, f = trimesh.proximity.closest_point(mesh2d, points2d)
    c = c[d < 1e-8]
    f = f[d < 1e-8]
    # Homogeneous coordinates
    c[:, 2] = 1
    p = []
    for ci, fi in zip(c, f):
        R = np.ones((3, 3))
        R[0] = u[mesh.faces[fi]]
        R[1] = v[mesh.faces[fi]]
        bary = np.linalg.solve(R, ci)
        p.append(mesh.vertices[mesh.faces[fi], :].T @ bary)

    return np.array(p)


#%% Determine 2D parameterization and plot coordinate function on the 3D mesh
from mayavi import mlab
from bfieldtools.viz import plot_data_on_vertices, plot_mesh, plot_data_on_faces

mesh = load_example_mesh("meg_helmet", process=False)
u, v, mesh2d = flatten_mesh(mesh, _lambda=1.0)

plot_data_on_vertices(mesh, u, ncolors=15)
plot_data_on_vertices(mesh, v, ncolors=15)

#%% Plot flattened mesh and area distortion on faces
plot_data_on_faces(mesh2d, mesh2d.area_faces / mesh.area_faces)

#%% Plot gradient of the two coordinate functions and the cosine of the angle between the gradients
from bfieldtools.mesh_calculus import gradient

gx = gradient(u, mesh)
gy = gradient(v, mesh)
cos = np.sum(gx * gy, axis=0) / (
    np.linalg.norm(gx, axis=0) * np.linalg.norm(gy, axis=0)
)
plot_data_on_faces(mesh, cos, vmin=-1, vmax=1)
mlab.quiver3d(*mesh.triangles_center.T, *gx, color=(1, 0, 0), mode="arrow")
mlab.quiver3d(*mesh.triangles_center.T, *gy, color=(0, 0, 1), mode="arrow")


#%% Map hexagonal grid from 2d to the 3D mesh
d = np.sqrt(3 / 4)
m = np.min((u.min(), v.min()))
mm = np.min((u.max(), v.max()))
xx = np.linspace(m * 1.05, mm * 1.05, 12)
yy = np.linspace(m * 1.05, mm * 1.05, 12) * d
p = np.array(np.meshgrid(xx, yy, 0, indexing="ij"))
p[0, :, ::2] += (xx[1] - xx[0]) * d / 2

p = p.reshape(3, -1).T

pp = plane2mesh(p, mesh, u, v)

plot_data_on_vertices(mesh, u, ncolors=15)
mlab.points3d(*pp.T, scale_factor=0.01)
