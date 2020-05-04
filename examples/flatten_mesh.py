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


mesh = load_example_mesh("meg_helmet", process=False)

N = mesh.vertices.shape[0]

L = laplacian_matrix(mesh)
M = mass_matrix(mesh)
Ac = A_matrix_complex(mesh)
vals, uv = eigsh(-L.T - Ac, 6, M, which="LM", sigma=0)

# Coordinates with inital phase
x = uv[:, 1].real
y = uv[:, 1].imag

#%% Determine "phase" by matching the coordinate function with mesh coordinates
theta = np.linspace(0, 2 * np.pi, 50)
xx = np.real(np.exp(1j * theta)[:, None] * uv[:, 1])
# plt.plot(np.sum(mesh.vertices[:,0]*xx, axis=1))
ii = np.argmax(np.sum(mesh.vertices[:, 0] * xx, axis=1))

theta = theta[ii]
x = np.real(np.exp(1j * theta) * uv[:, 1])
y = np.imag(np.exp(1j * theta) * uv[:, 1])
#%%
from mayavi import mlab
from bfieldtools.viz import plot_data_on_vertices, plot_mesh

plot_data_on_vertices(mesh, x, ncolors=15)
plot_data_on_vertices(mesh, y, ncolors=15)

#%%
mlab.figure()
mlab.triangular_mesh(
    x,
    y,
    np.zeros_like(x),
    mesh.faces,
    scalars=x,
    representation="wireframe",
    colormap="bwr",
)

#%% Plot gradient of the two coordinate functions
from bfieldtools.mesh_calculus import gradient

gx = gradient(x, mesh)
gy = gradient(y, mesh)
plot_mesh(mesh)
mlab.quiver3d(*mesh.triangles_center.T, *gx, color=(1, 0, 0), mode="arrow")
mlab.quiver3d(*mesh.triangles_center.T, *gy, color=(0, 0, 1), mode="arrow")
