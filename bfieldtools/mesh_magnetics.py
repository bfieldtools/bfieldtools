'''
Contains functions for calculating the coupling of surface current density in a
triangle mesh to magnetic field as well as scalar and vector potentials.

'''

import numpy as np
#from numba import jit
import time

from .utils import get_quad_points
from .mesh_calculus import gradient_matrix
from .integrals import triangle_potential_dipole_linear, triangle_potential_uniform



def magnetic_field_coupling(mesh, r, Nchunks=None, quad_degree=1):
    '''
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object describing mesh
    r: target points (Np, 3)
    quad_degree: int >= 1
        Quadrature degree (Dunavant scheme) to use.

    Returns
    -------
    C: (Np, 3, Nvertices) array
        Coupling matrix for surface current in the mesh to the evaluation points)

    '''
    mu0 = 4 * np.pi * 1e-7
    coef = mu0 / (4 * np.pi)

    print('Computing magnetic field coupling matrix, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    w_quad, r_quad = get_quad_points(mesh.vertices, mesh.faces, method='dunavant_0'+str(quad_degree))

    # Rotated gradients (currents)
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)

    # Initialize C-matrix
    n_target_points = len(r)
    n_verts = len(mesh.vertices)
    C = np.zeros((n_target_points, n_verts, 3))

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

    for n in range(Nchunks):
        # Diffence vectors (Neval, Ntri, Nquad, 3)
        RR = r_quad[None, :, :, :] - r[n::Nchunks, None, None, :]

        # RR/norm(RR)**3 "Gradient of Green's function"
        g = - RR/((np.linalg.norm(RR, axis=-1)**3)[:, :, :, None])

        # Sum over quad points and multiply by triangle area
        g = (g*w_quad[:, None]).sum(axis=-2)
        g *= mesh.area_faces[:, None]

        # Cross product RR/norm(RR)
        C[n::Nchunks, :, 0] = g[:, :, 2] @ Gy - g[:, :, 1] @ Gz
        C[n::Nchunks, :, 1] = g[:, :, 0] @ Gz - g[:, :, 2] @ Gx
        C[n::Nchunks, :, 2] = g[:, :, 1] @ Gx - g[:, :, 0] @ Gy


    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    return coef * np.moveaxis(C, 2, 1)


def magnetic_field_coupling_analytic(mesh, r, Nchunks=None):
    '''
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object describing the mesh
    r: target points (Np, 3)

    Returns
    -------
    C: (Np, 3, Nvertices) array
        Coupling matrix for surface current in the mesh to the evaluation points)

    '''
    from .integrals import omega, gamma0
    coef = 1e-7

    print('Computing magnetic field coupling matrix analytically, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

#    ta = mesh.area_faces
    tn = mesh.face_normals

    # Nfaces, 3, 3
    rfaces = mesh.vertices[mesh.faces]
    # Neval, Nfaces, xyz
    coeffs = np.zeros((r.shape[0:1] + mesh.faces.shape[:-1] + (3,)))
    # Nfaces, Nedges, 3
    edges = np.roll(rfaces, 1, -2) - np.roll(rfaces, 2, -2)
#    grad = np.cross(tn[:, None, :], edges, axis=-1)/(2*ta[:, None, None])
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=False)
    Rx, Ry, Rz = gradient_matrix(mesh, rotated=True)
    solid_angle = np.zeros((r.shape[0], rfaces.shape[0]))
    # Calculate potentials and related coefficients
    for n in range(Nchunks):
        RRchunk = r[n::Nchunks, None, None, :] - rfaces[None, :, :, :]
        # Neval, Nfaces, Nedges
#        result = -np.sum(np.sum(gamma0(RRchunk)[..., None]*edges,
#                               axis=-2)[...,None, :]*edges, axis=-1)
#        result *= 1/(2*ta[..., :, None])
        solid_angle[n::Nchunks] = omega(RRchunk)
        coeffs[n::Nchunks] = -np.einsum('...i,...ik->...k', gamma0(RRchunk), edges)

#    # Accumulate the elements
    C = np.zeros((3, r.shape[0], mesh.vertices.shape[0]))
#    for ind_f, f in enumerate(mesh.faces):
#        C[:, f, :] += coeffs[:, ind_f, :, None]*tn[ind_f]
#        C[:, f, :] += -solid_angle[:, ind_f:ind_f+1, None]*grad[ind_f]

    G = (Gx, Gy, Gz)
    for i in range(3):
        cc = coeffs*tn[:, i, None]
        C[i] += cc[:,:,0] @ Rx + cc[:,:,1] @ Ry + cc[:,:,2] @ Rz
        C[i] += -solid_angle @ G[i]

    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    return coef * np.moveaxis(C, 0, 1)


def scalar_potential_coupling(mesh, r, Nchunks=None):
    """
    Compute scalar potential matrix from linear stream functions
    using analytic integral

    Parameters
    ----------

    mesh: Trimesh mesh object describing mesh
    r: target points (Np, 3)

    Returns
    -------
    U: (Np, 3, Nvertices) array
        Coupling matrix for surface current in the mesh to the evaluation points
    """

    coeff = 1e-7  # mu_0/(4*pi)

    print('Computing scalar potential coupling matrix, %d vertices by %d target points... '%(len(mesh.vertices), len(r)), end='')
    start = time.time()

    # Source and eval locations
    R1 = mesh.vertices[mesh.faces]
    R2 = r

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

    R2chunks = np.array_split(R2, Nchunks, axis=0)
    i0=0
    Uf = np.zeros((R2.shape[0], mesh.faces.shape[0], 3))
    for R2chunk in R2chunks:
        RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
        i1 = i0+RRchunk.shape[0]
        Pi = triangle_potential_dipole_linear(RRchunk, mesh.face_normals,
                                             mesh.area_faces)
        Uf[i0:i1] = Pi
        i0=i1

#     Accumulate the elements
    Uv = np.zeros((R2.shape[0], mesh.vertices.shape[0]))
    for ind_f, f in enumerate(mesh.faces):
        Uv[:, f] += Uf[:, ind_f]


    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    return Uv*coeff


def vector_potential_coupling(mesh, r, Nchunks=None):
    """
    Compute vector potential matrices (one for each coordinate)
    from linear stream functions using analytic integral
    Parameters
    ----------

    mesh: Trimesh mesh object describing mesh
    r: target points (Np, 3)

    Returns
    -------
    A: (Np, 3, Nvertices) array
        Coupling matrix for surface current in the mesh to the evaluation points
    """

    coeff = 1e-7  # mu_0/(4*pi)

    # Source and eval locations
    R1 = mesh.vertices[mesh.faces]
    R2 = r

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0]//100
        else:
            Nchunks = 1

    R2chunks = np.array_split(R2, Nchunks, axis=0)
    i0=0
    Af = np.zeros((R2.shape[0], mesh.faces.shape[0]))
    print('Computing potential matrix')

    for R2chunk in R2chunks:
        RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
        i1 = i0+RRchunk.shape[0]
        Pi = triangle_potential_uniform(RRchunk, mesh.face_normals, False)
        Af[i0:i1] = Pi
#        print((100*i1)//R2.shape[0], '% computed')
        i0=i1

    # Rotated gradients (currents)
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
    # Accumulate the elements
    Av = np.array([Af @ Gx, Af @ Gy, Af @ Gz])

    return Av*coeff