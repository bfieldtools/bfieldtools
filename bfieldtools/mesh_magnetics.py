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
from .integrals import triangle_potential_approx, potential_dipoles



def magnetic_field_coupling(mesh, r, Nchunks=None, quad_degree=1, analytic=False):
    '''
    Given 'mesh', computes the "C matrix" which gives the magnetic field at
    target points 'r' due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evalution points
    quad_degree: int >= 1
        Quadrature degree (Dunavant scheme) to use.
    analytic: boolean
        compute field using analytic formula (True) or quadrature (False)

    Returns
    -------
    C: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to B-field at the evaluation points

    '''

    if analytic:
        return magnetic_field_coupling_analytic(mesh, r, Nchunks)

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

    C *= coef
#    return np.moveaxis(C, 2, 1)
    return np.swapaxes(C, 2, 1)



def magnetic_field_coupling_analytic_old(mesh, r, Nchunks=None):
    '''
    DEPRECATED
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object
            mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points

    Returns
    -------
    C: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to B-field at the evaluation points

    DEPRECATED
    '''
    from .integrals_old import omega, gamma0
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

    C *= coef
#    return np.moveaxis(C, 0, 1)
    return np.swapaxes(C, 0, 1)



def magnetic_field_coupling_analytic(mesh, r, Nchunks=None):
    '''
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh
    using analytic formulas.

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points
    Nchunks: int
        number of chunks used in the calculation for saving memory

    Returns
    -------
    C: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to B-field at the evaluation points

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

    ta = mesh.area_faces
    tn = mesh.face_normals

    # Nfaces, 3, 3
    rfaces = mesh.vertices[mesh.faces]

    # Calculate potentials and related coefficients
    gamma_terms = np.zeros((r.shape[0], mesh.faces.shape[0], 3))
    omega_terms = np.zeros((r.shape[0], mesh.faces.shape[0]))
    # Edges Nfaces, 3, 3
    edges = np.roll(rfaces, 1, -2) - np.roll(rfaces, 2, -2)
    for n in range(Nchunks):
        RRchunk = r[n::Nchunks, None, None, :] - rfaces[None, :, :, :]
        # Neval, Nfaces, xyz
        gamma_terms[n::Nchunks] = -np.einsum('nfe,fei->nfi', gamma0(RRchunk), edges)
        omega_terms[n::Nchunks] = omega(RRchunk)

    # 3 (Nfaces, Nverts) sparse matrices
    G = gradient_matrix(mesh, rotated=False)
    R = gradient_matrix(mesh, rotated=True)
    C = np.zeros((3, r.shape[0], mesh.vertices.shape[0]))

    # Accumulate elements by sparse matrix products for x,y, and z components
    for fcomp in range(3):
        C[fcomp] = omega_terms @ G[fcomp]
        # Accumulate gamma terms for each vertex in the triangle
        for gcomp in range(3):
            # Edge @ Rotated_gradient "==" c_coeff
            # Multiplying with R-matrices takes care of c_coeff calculation
            # and accumulation to right vertex
            C[fcomp] += (tn[:, fcomp]*gamma_terms[:, :, gcomp]) @ R[gcomp]

    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    C *= coef
#    return np.moveaxis(C, 0, 1)
    return np.swapaxes(C, 0, 1)



def scalar_potential_coupling(mesh, r, Nchunks=None, multiply_coeff=True,
                              approx_far=False, margin=3):
    """
    Coupling matrix corresponding to a mapping from a stream function
    to scalar potential using analytic integrals

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points
    Nchunks: int
        number of chunks used in the calculation for saving memory
    multiply_coeff: boolean
        If True, multiply result by mu_0/(4*pi)
    approx_far: boolean,
        approximate the potential using simple quadrature
        (see integrals.potential_dipoles) for points far from the source triangles
    margin: float
        cut-off distance for "far" points measured in mean triangle side length.

    Returns
    -------
    U: ndarray (Np, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to scalar potential at the evaluation points
    """

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
        if approx_far:
            near, far = _split_by_distance(mesh, RRchunk, margin)
            Uf[i0:i1, near] = triangle_potential_dipole_linear(RRchunk[:, near],
                                                               mesh.face_normals[near],
                                                               mesh.area_faces[near])
            Uf[i0:i1, far] = potential_dipoles(RRchunk[:, far],
                                               mesh.face_normals[far],
                                               mesh.area_faces[far])
        else:
            Uf[i0:i1] = triangle_potential_dipole_linear(RRchunk, mesh.face_normals,
                                                 mesh.area_faces)
        i0=i1

#     Accumulate the elements
    Uv = np.zeros((R2.shape[0], mesh.vertices.shape[0]))
    for ind_f, f in enumerate(mesh.faces):
        Uv[:, f] += Uf[:, ind_f]


    duration = time.time() - start
    print('took %.2f seconds.'%duration)

    if multiply_coeff:
        coeff = 1e-7  # mu_0/(4*pi)
    else:
        coeff = 1
    return Uv*coeff


def vector_potential_coupling(mesh, r, Nchunks=None, approx_far=True, margin=2):
    """
    Compute vector potential coupling matrices
    from linear stream functions using analytic integral

    Parameters
    ----------

    mesh: Trimesh mesh object
        mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points
    approx_far: Boolean (True)
        If True, use approximate calculation for triangles that
        far from the source triangles using a simple quadrature
        (see integrals.triangle_potential_approx)
    margin: float,
        cut-off distance for "far" points measured in mean triangle side length.

    Returns
    -------
    A: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to vector potential at the evaluation points
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
    print('Computing 1/r-potential matrix')

    for chunk_idx, R2chunk in enumerate(R2chunks):
#        print('Computing chunk %d/%d'%(chunk_idx+1, Nchunks))
        RRchunk = R2chunk[:, None, None, :] - R1[None, :, :, :]
        i1 = i0+RRchunk.shape[0]
        if approx_far:
            near, far = _split_by_distance(mesh, RRchunk, margin)
            Af[i0:i1, near] = triangle_potential_uniform(RRchunk[:, near], mesh.face_normals[near], False)
            Af[i0:i1, far] = triangle_potential_approx(RRchunk[:, far], mesh.area_faces[far], reg=0)
        else:
            Af[i0:i1] = triangle_potential_uniform(RRchunk, mesh.face_normals, False)
#        print((100*i1)//R2.shape[0], '% computed')
        i0=i1

    #Free some memory by deleting old variables
    del R1, R2, RRchunk, R2chunks

    # Rotated gradients (currents)
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
    # Accumulate the elements
    Av = np.array([Af @ Gx, Af @ Gy, Af @ Gz])

    return Av*coeff


def _split_by_distance(mesh, RR, margin=3):
    avg_sidelength = 8/np.sqrt(3)*np.mean(mesh.area_faces[::100])
#    np.mean(np.linalg.norm(np.diff(mesh.vertices[mesh.edges[::1000]], axis=1), axis=-1))

    RRnorm = np.linalg.norm(RR, axis=-1)
    near = np.nonzero(np.min(RRnorm, axis=(0, 2)) < avg_sidelength * margin)[0]

    far = np.setdiff1d(np.arange(0, len(mesh.faces)), near, assume_unique=True)

#    print('near: %d, far: %d'%(len(near), len(far)))
    return near, far


