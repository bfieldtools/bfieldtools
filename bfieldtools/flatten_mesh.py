"""
functions for flattening a mesh with a boundary
"""

__all__ = [
    "A_matrix_complex",
    "eigen_complex_laplacian",
    "flatten_mesh",
    "mesh2plane",
    "plane2mesh",
]

import numpy as np
from bfieldtools.mesh_calculus import laplacian_matrix, mass_matrix
from scipy.sparse.linalg import eigsh
from scipy.sparse import csc_matrix
import trimesh


def A_matrix_complex(mesh):
    """
    Area matrix for complex parametrization
    See this course:
    https://www.cs.cmu.edu/~kmcrane/Projects/DDG/
    section 7.4.
    
    CURRENTLY ONLY FOR A MESH WITH A SINGLE BOUNDARY (no holes)

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


def eigen_complex_laplacian(mesh, Nc, _lambda):
    L = laplacian_matrix(mesh)
    M = mass_matrix(mesh)
    Ac = A_matrix_complex(mesh)
    vals, uv = eigsh(-0.5 * L.T - _lambda * Ac, Nc, M, which="LM", sigma=0)
    return vals, uv


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

    _lambda <= 1.0
    _lambda == 1.0 => conformal mapping
    _lambda == 0.5 =>  not conformal but less area distortion
    _lambda --> 0 mapping becomes denegerate (real==imag)
    _lambda > 1 (e.g. 1.01-1.1) folding effects
    """
    vals, uv = eigen_complex_laplacian(mesh, 2, _lambda)

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
    """
    Map points from 3D to 2D on the u,v -plane

    Parameters
    ----------
    points3d : nparray
        points close to the  mesh (N, 3)
    mesh : Trimesh object
        3D mesh with the u,v -parametrisation
    u : ndarray
        u coordinate (N,)
    v : ndarray
        v coordinate (N,)

    Returns
    -------
    ndarray
       uv-coordinates of the 3D points

    """
    c, d, f = trimesh.proximity.closest_point(mesh, points3d)
    tris = mesh.vertices[mesh.faces[f]]
    barys = trimesh.triangles.points_to_barycentric(tris, c)
    print(barys)
    p1 = np.sum(u[mesh.faces[f]] * barys, axis=1)
    p2 = np.sum(v[mesh.faces[f]] * barys, axis=1)
    return np.array([p1, p2]).T


def plane2mesh(points2d, mesh, u, v):
    """
    Map point from 2D to 3D on the mesh

    Parameters
    ----------
    points3d : nparray
        points close to the  mesh (N, 3)
    mesh : Trimesh object
        3D mesh with the u,v -parametrisation
    u : ndarray
        u coordinate (N,)
    v : ndarray
        v coordinate (N,)

    Returns
    -------
    ndarray
       uv-coordinates of the 3D points
    """
    mesh2d = trimesh.Trimesh(np.array([u, v, 0 * u]).T, mesh.faces)
    c, d, f = trimesh.proximity.closest_point(mesh2d, points2d)
    # Leave points outside the mesh area
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
