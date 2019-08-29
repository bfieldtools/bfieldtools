'''
Contains functions for computing the inductance matrices of triangle surface meshes,
including both self- and mutual-inductance.
'''

import numpy as np
from .utils import (get_quad_points, assemble_matrix_chunk, assemble_matrix2)


def gamma0(R, reg=1e-13, symmetrize=True):
    """ Integrals over the edges of a triangle called gamma_0 (line charge potentials).

        **NOTE: MAY NOT BE VERY PRECISE FOR POINTS DIRECTLY AT TRIANGLE
        EDGES.**

        Parameters
        ----------
        R : (N, 3, 3) array of points (Neval, Nverts, xyz)

        Returns
        -------
        res: array (Neval, Nverts)
            The analytic integrals for each vertex/edge



    """
    diffs = np.roll(R, 1, -2) - np.roll(R, 2, -2)
    dotprods1 = np.sum(np.roll(R, 1, -2)*diffs, axis=-1)
    dotprods2 = np.sum(np.roll(R, 2, -2)*diffs, axis=-1)
    dn = np.linalg.norm(diffs, axis=-1)
    del diffs
    n = np.linalg.norm(R, axis=-1)
    # Regularize s.t. neither the denominator or the numerator can be zero
    # Avoid numerical issues directly at the edge
    res = np.log((np.roll(n, 2, -1)*dn + dotprods2 + reg)
                 / (np.roll(n, 1, -1)*dn + dotprods1 + reg))

    # Symmetrize the result since on the negative extension of the edge
    # there's division of two small values resulting numerical instabilities
    # (also incompatible with adding the reg value)
    if symmetrize:
        res2 = -np.log((np.roll(n, 2, -1)*dn - dotprods2 + reg)
                       / (np.roll(n, 1, -1)*dn - dotprods1 + reg))
        res = np.where(dotprods1+dotprods2 > 0, res, res2)
    res /= dn
    return res


def omega(R):
    """ Calculate the solid angle of a triangles using
        A. Van Oosterom and J. Strackee
        IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING,
        VOL. BME-30, NO. 2, 1983

        Parameters
        ----------
        R : (N, ..., 3, 3) array of points (Neval, ..., Nverts, xyz)
            Points correspond to relative coordinates (x,y,z) of
            N triangles/evaluation points for
            the 3 corners of the triangles/triangle.

            Neval can be number of evaluation points for the same triangle
            or number of triangles for the same evaluation points

        Returns
        -------
        sa: (Neval, 3) ???
            Solid angles of triangle(s) at evaluation points
    """
    # Distances
    d = np.linalg.norm(R, axis=-1)
    # Scalar triple products
    stp = np.linalg.det(R)
    # Denominator
    denom = np.prod(d, axis=-1)
    for i in range(3):
        j = (i+1) % 3
        k = (i+2) % 3
        denom += np.sum(R[..., i, :]*R[..., j, :], axis=-1)*d[..., k]
    # Solid angles
    sa = 2*np.arctan2(stp, denom)
    return sa


def triangle_potential(R, tn, planar=False):
    """ 1/r potential of a uniform triangle

        Parameters
        ----------

        R : (N, (Ntri), 3, 3) array
            Displacement vectors (Neval, ...., Ntri_verts, xyz)
        tn : ((Ntri), 3) array
            Triangle normals (Ntri, dir)
        planar: boolean
            If True, use planar geometry assumption for speed

        Returns
        -------
        result:  array
            Resultant 1/r potential

    """
    if len(R.shape) > 3:
        tn_ax = tn[:, None, :]
    else:
        tn_ax = tn
    summands = gamma0(R)*np.sum(tn_ax*np.cross(np.roll(R, 1, -2),
                                               np.roll(R, 2, -2), axis=-1), axis=-1)
    if not planar:
        csigned = np.sum(np.take(R, 0, -2)*tn, axis=-1)
        result = np.sum(summands, axis=-1) - csigned*omega(R)
    else:
        print('Assuming all the triangles are in the same plane!')
        result = np.sum(summands, axis=-1)
    return result


def triangle_potential_approx(R, ta, planar=False, reg=1e-12):
    """ 1/r potential of a uniform triangle using centroid approximation

        Calculates 1/R potentials for triangle centroids
        (The singular at the centroid is handled with the very small
        reg value, but anyway the values close to the centroid are inexact)

        Parameters
        ----------
        R : (N, (Ntri), 3, 3) array
            Displacement vectors (Neval, ...., Ntri_verts, xyz)
        ta : (Ntri) array
            Triangle areas
        planar: boolean
            If True, use planar geometry assumption for speed
        reg: float
            Regularization value used in approximation

        Returns
        -------
        result: array
            Resultant 1/r potential

    """
    result = 1/(np.linalg.norm(np.mean(R, axis=-2), axis=-1) + reg)*ta
    return result

def mutual_inductance_matrix(mesh1, mesh2, planar=False):
    """ Calculate a mutual inductance matrix for hat basis functions
        (stream functions) between two surface meshes

        Parameters
        ----------

        mesh1: Trimesh mesh object for mesh 1
        mesh2: Trimesh mesh object for mesh 2
        planar: boolean
            If True, use planar assumption when calculating

        Returns
        -------
        M: (Nvertices1 x Nvertices2) array
            Mutual inductance matrix between mesh1 and mesh2

    """
    R = mesh1.vertices[mesh1.faces]  # Nt x 3 (corners) x 3 (xyz)
    # Calculate quadrature points
    weights, quadpoints = get_quad_points(mesh2.vertices, mesh2.faces, 'centroid')
    # Nt x Nquad x  3 (x,y,z)

    RR = quadpoints[:, :, None, None, :] - R[None, None, :, :, :]
    print('Calculating potentials')
    pots = triangle_potential(RR, mesh1.face_normals, planar=planar) # Ntri_eval, Nquad, Ntri_source
    pots = np.sum(pots * weights[None, :, None], axis=1) # Ntri_eval, Ntri_source

    # Calculate edge vectors for each triangle
    edges1 = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)
    edges2 = np.roll(mesh2.vertices[mesh2.faces], 1, -2) - np.roll(mesh2.vertices[mesh2.faces], 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)

    tri_data = np.sum(edges1[None, :, None, :, :] * edges2[:, None, :, None, :], axis=-1) # i,j,k,l
    tri_data /= (mesh2.area_faces[:, None] * mesh1.area_faces[None, :] * 4)[:, :, None, None]
    tri_data *= (mesh2.area_faces[:, None] * pots)[:, :, None, None]
    print('Inserting stuff into M-matrix')

    M = assemble_matrix2(mesh1.faces, mesh2.faces, mesh1.vertices.shape[0], mesh2.vertices.shape[0], tri_data)
    return M * 1e-7


def self_inductance_matrix(mesh, planar=False, Nchunks=1, approx=False):
    """ Calculate a self inductance matrix for hat basis functions
        (stream functions) in the triangular mesh described by

        Parameters
        ----------
        mesh: Trimesh mesh object

        Returns
        -------
        M: (Nvertices x Nvertices) array
            Self.inductance matrix of `mesh`
    """
    R = mesh.vertices[mesh.faces]  # Nt x 3 (corners) x 3 (xyz)
    # Calculate edge vectors for each triangle
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2)  # Nt x 3 (edges) x 3 (x,y,z)
    # Calculate quadrature points
    weights, quadpoints = get_quad_points(mesh.vertices, mesh.faces, 'centroid')
    # Nt x Nquad x  3 (x,y,z)


    M = np.zeros((mesh.vertices.shape[0], mesh.vertices.shape[0]))

    for n in range(Nchunks):
        RR = quadpoints[n::Nchunks, :, None, None, :] - R[None, None, :, :, :]

        # Loop evaluation triangles (quadpoints) in chunks
        # Init not needed
#        tri_data = np.zeros((len(quadpoints[n::Nchunks]), edges.shape[0],
#                             edges.shape[1], edges.shape[1]))


        print('Calculating potentials, chunk %d/%d' % (n + 1, Nchunks))
        if approx:
            pots = triangle_potential_approx(RR, mesh.area_faces) # Ntri_eval, Nquad, Ntri_source
            # Recalculate the nearby potentials with analytical formula
            # "diagonal" indices
#           i0 = i1 = np.arange(RR.shape[0])
            # Find 1-neighbourhood around vertices
            fsparse = mesh.faces_sparse.tocsc()[:,n::Nchunks]
            nb_inds = np.nonzero((fsparse.T @ fsparse).toarray())
            i0=nb_inds[0]
            i1=nb_inds[1]
            # RR  and normals of the neighbourhoods
            RR_nb = RR[:,:,n::Nchunks,:,:][i0,:,i1,:,:]
            RR_nb = np.moveaxis(RR_nb, 0, 1)
            n_nb = mesh.face_normals[n::Nchunks][i1]
            pots[:, :, n::Nchunks][i0, :, i1] = triangle_potential(RR_nb, n_nb, planar=planar).T # Ntri_eval, Nquad, Ntri_source
        else:
            pots = triangle_potential(RR, mesh.face_normals, planar=planar) # Ntri_eval, Nquad, Ntri_source
        pots = np.sum(pots*weights[None, :, None], axis=1) # Ntri_eval, Ntri_source

        tri_data = np.sum(edges[None, :, None, :, :] * edges[n::Nchunks, None, :, None, :], axis=-1) # i,j,k,l
        tri_data /= (mesh.area_faces[n::Nchunks, None] * mesh.area_faces[None, :] * 4)[:, :, None, None]
        tri_data *= (mesh.area_faces[n::Nchunks, None] * pots)[:, :, None, None]

        M += assemble_matrix_chunk(mesh.faces, mesh.vertices.shape[0], tri_data, n, Nchunks)

    return M * 1e-7


