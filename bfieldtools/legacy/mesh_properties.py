from ..mesh_impedance import mutual_inductance_matrix

import numpy as np


def self_inductance_matrix(
    mesh, Nchunks=None, quad_degree=2, approx_far=True, margin=2, chunk_clusters=False,
):
    """ Calculate a self inductance matrix for hat basis functions
        (stream functions) in the triangular mesh described by

        Parameters
        ----------
        mesh: Trimesh mesh object
        Nchunks: int
            Number of serial chunks to divide the computation into
        quad_degree: int >= 1
            Quadrature degree (Dunavant scheme) to use. Self-inductance requires higher degree than mutual inductance
        approx_far: Boolean (True)
            If True, use approximate calculation for triangles that
            far from the source triangles using a simple quadrature
            (see integrals.triangle_potential_approx)
        margin: float
            Cut-off distance for "far" points measured in mean triangle side length

        Returns
        -------
        M: (Nvertices x Nvertices) array
            Self.inductance matrix of `mesh`
    """
    if quad_degree <= 2:
        print(
            "Computing self-inductance matrix using rough quadrature (degree=%d).\
              For higher accuracy, set quad_degree to 4 or more."
            % quad_degree
        )

    return mutual_inductance_matrix(
        mesh,
        mesh,
        Nchunks=Nchunks,
        quad_degree=quad_degree,
        approx_far=approx_far,
        margin=margin,
    )


def triangle_self_coupling(mesh):
    """
    Self-coupling integrated analytically. Implemented based on
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=475946

    Parameters
    ----------
    mesh: Trimesh mesh object

    Returns
    -------
    self_coupling: array (N_triangles, )
        triangle self-coupling
    """

    tri_points = mesh.vertices[mesh.faces]

    r1 = tri_points[:, 0, :]
    r2 = tri_points[:, 1, :]
    r3 = tri_points[:, 2, :]

    a = np.einsum("ij,ij->i", r3 - r1, r3 - r1)
    b = np.einsum("ij,ij->i", r3 - r1, r3 - r2)
    c = np.einsum("ij,ij->i", r3 - r2, r3 - r2)

    sa = np.sqrt(a)
    sc = np.sqrt(c)
    ss = np.sqrt(a - 2 * b + c)
    sac = np.sqrt(a * c)

    self_coupling = (1 * (4 * mesh.area_faces ** 2)) * (
        1
        / (6 * sa)
        * np.log(((a - b + sa * ss) * (b + sac)) / ((-b + sac) * (-a + b + sa * ss)))
        + 1
        / (6 * sc)
        * np.log(((b + sac) * (-b + c + sc * ss)) / ((b - c + sc * ss) * (-b + sac)))
        + 1
        / (6 * ss)
        * np.log(
            ((a - b + sa * ss) * (-b + c + sc * ss))
            / ((b - c + sc * ss) * (-a + b + sa * ss))
        )
    )

    return self_coupling
