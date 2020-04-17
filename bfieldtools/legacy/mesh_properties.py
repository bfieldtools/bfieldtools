from ..mesh_properties import mutual_inductance_matrix


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
