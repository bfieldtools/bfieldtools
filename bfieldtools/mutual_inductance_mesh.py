'''
Contains functions for computing the inductance matrices of triangle surface meshes,
including both self- and mutual-inductance.
'''

import numpy as np
from .utils import (get_quad_points, assemble_matrix_chunk, assemble_matrix2)
from .integrals import triangle_potential_uniform as triangle_potential
from .integrals import triangle_potential_approx



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


