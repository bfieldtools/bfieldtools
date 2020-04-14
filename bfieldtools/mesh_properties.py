"""
Contains functions for computing the inductance matrices of triangle surface meshes,
including both self- and mutual-inductance.
"""
from psutil import virtual_memory
import numpy as np

from .utils import get_quad_points, get_line_quad_points
from .mesh_magnetics import vector_potential_coupling
from .mesh_calculus import gradient_matrix, laplacian_matrix


def resistance_matrix(mesh, sheet_resistance):
    """ Resistance matrix

        Parameters
        ----------
        mesh: Trimesh mesh object
        sheet_resistance: (N_faces) array or scalar
                            "1/sigma*d" constant resistance for each face (or all faces if scalar)
        Returns
        -------
        R: (Nvertices x Nvertices) array
            resistance matrix of `mesh`
    """
    return -laplacian_matrix(mesh, sheet_resistance)


def self_inductance_matrix(
    mesh, Nchunks=None, quad_degree=2, approx_far=True, margin=2, chunk_clusters=False
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


def mutual_inductance_matrix(
    mesh1, mesh2, Nchunks=None, quad_degree=1, approx_far=True, margin=2
):
    """ Calculate a mutual inductance matrix for hat basis functions
        (stream functions) between two surface meshes

        Parameters
        ----------

        mesh1: Trimesh mesh object for mesh 1
        mesh2: Trimesh mesh object for mesh 2
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
        M: (Nvertices1 x Nvertices2) array
            Mutual inductance matrix between mesh1 and mesh2

    """

    if Nchunks is None:
        # Available RAM in megabytes
        mem = virtual_memory().available >> 20

        # Estimate of memory usage in megabytes for a single chunk, when quad_degree=2 (very close with quad_degree=1)
        mem_use = 0.033 * (len(mesh1.vertices) * len(mesh2.vertices)) ** 0.86

        print(
            "Estimating %d MiB required for %d by %d vertices..."
            % (mem_use, len(mesh1.vertices), len(mesh2.vertices))
        )

        # Chunk computation so that available memory is sufficient
        Nchunks = int(np.ceil(mem_use / mem))

        if approx_far:
            Nchunks *= 20
            print(
                "Computing inductance matrix in %d chunks (%d MiB memory free),\
                  when approx_far=True using more chunks is faster..."
                % (Nchunks, mem)
            )
        else:
            print(
                "Computing inductance matrix in %d chunks since %d MiB memory is available..."
                % (Nchunks, mem)
            )

    # Calculate quadrature points
    weights, quadpoints = get_quad_points(
        mesh2.vertices, mesh2.faces, "dunavant_0" + str(quad_degree)
    )
    # Nt x Nquad x  3 (x,y,z)

    # Compute vector potential to quadrature points
    Nw = len(weights)
    Nt = len(mesh2.faces)
    Nv = len(mesh1.vertices)

    A = vector_potential_coupling(
        mesh1,
        quadpoints.reshape(-1, 3),
        Nchunks=Nchunks,
        approx_far=approx_far,
        margin=margin,
    ).reshape(3, Nt, Nw, Nv)

    # Integrate over the triangles (current patterns are constant over triangles)
    A = np.sum(A * weights[None, None, :, None], axis=2)
    A *= mesh2.area_faces[None, :, None]

    # Dot product with current patterns and sum over triangle neighbourhoods
    Gx, Gy, Gz = gradient_matrix(mesh2, rotated=True)
    M = A[0].T @ Gx + A[1].T @ Gy + A[2].T @ Gz

    return M


def triangle_self_coupling(mesh):
    """
    Self-coupling integrated analytically. Implemented based on
    Poole, M.S., 2007. Improved equipment and techniques for dynamic shimming in high field MRI (Doctoral dissertation, University of Nottingham.). page 72.


    Self-coupling can be integrated analytically using different notation
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=475946
    """

    tri_points = mesh.vertices[mesh.faces]

    r1 = tri_points[:, 0, :]
    r2 = tri_points[:, 1, :]
    r3 = tri_points[:, 2, :]

    a = np.dot(r3 - r1, r3 - r1)
    b = np.dot(r3 - r1, r3 - r2)
    c = np.dot(r3 - r2, r3 - r2)

    a = np.einsum("ij,ij->i", r3 - r1, r3 - r1)
    b = np.einsum("ij,ij->i", r3 - r1, r3 - r2)
    c = np.einsum("ij,ij->i", r3 - r2, r3 - r2)

    sa = np.sqrt(a)
    sc = np.sqrt(c)
    ss = np.sqrt(a - 2 * b + c)
    sac = np.sqrt(a * c)

    integral = (1 * (4 * mesh.area_faces ** 2)) * (
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

    return integral


def mesh2line_mutual_inductance(mesh, line_vertices, quad_degree=3):
    """
    Mutual inductance of a closed line segment loop (last segment connecting to first)
    and a triangle mesh

    Parameters
    ----------
    mesh: Trimesh mesh object
    line_vertices: points connected in index order (N_points, 3)

    Returns
    -------

    M: mutual inductance vector with shape (N_vertices,)


    """

    # Calculate quadrature points
    weights, quadpoints = get_line_quad_points(
        line_vertices, "gauss_legendre", quad_degree
    )
    # Ne x Nquad x  3 (x,y,z)

    segments = np.roll(line_vertices, shift=-1, axis=0) - line_vertices

    # Compute vector potential to quadrature points
    Nw = len(weights)
    Nt = len(line_vertices)
    Nv = len(mesh.vertices)
    M = vector_potential_coupling(mesh, quadpoints.reshape(-1, 3)).reshape(
        3, Nt, Nw, Nv
    )

    # Integrate over quadrature points
    M = np.sum(M * weights[None, None, :, None], axis=2)

    # Scale by segment lengths, integrate over xyz-axis and segments
    M = np.sum(segments.T[:, :, None] * M, axis=(0, 1))

    return M
