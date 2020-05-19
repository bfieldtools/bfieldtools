"""

Functions for creating contours (isolines) of a scalar function defined on a triangle mesh surface. Also contains functions for modifying the generated contours.

"""

__all__ = ["scalar_contour", "simplify_contour"]

import numpy as np

from scipy.sparse import eye as speye
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve

from .mesh_calculus import gradient
from . import mesh_conductor


def scalar_contour(mesh, scalars, N_contours=10, contours=None, return_values=False):
    """
    Computes contour loops (isolines) for a scalar function defined on a mesh.
    The winding direction of the loops is defined according to the rotated gradient of the scalar

    Parameters
    ----------
    mesh: Trimesh object
        Trimesh object containing the mesh on which the scalar function is defined
    scalars: array-like or StreamFunction
        Vector containing the values of the scalar function at each of the mesh vertices.
        If StreamFunction, uses the vertex-wise values
    N_contours: int
        Number of contours to generate
    contours: array-like
        Optional argument for manual input of contour levels. Overrides `N_contours`
    return_values: Boolean
        If True, also return contour values

    Returns
    -------
    contour_polys: list
        list with length `N_contours`. Each list element is anumpy array containing the
        coordinats of each polygon vertex.
    contour_values: array-like
        Vector containing the scalar function value for each contour line,
        returned if return_values is True
    """

    if isinstance(scalars, mesh_conductor.StreamFunction):
        scalars = scalars.vert

    # Compute rotated gradient of scalars (e.g. current density from a stream function)
    g = gradient(scalars, mesh, rotated=True)

    if contours is None:
        # N evenly spaced contours that do not contain the point-like max and min contours
        contours = np.linspace(scalars.min(), scalars.max(), 2 * N_contours + 1)[1::2]

    edge_vals = scalars[mesh.edges_unique]
    contour_polys = []

    contour_values = []

    # Loop through
    for c in contours:

        # Get edges that contain the contour values
        edge_inds = (edge_vals.min(axis=1) <= c) * (edge_vals.max(axis=1) >= c)
        c_edge_vals = edge_vals[edge_inds]

        # Solve weights (barycentric coordinates) for each edge
        w0 = (c - c_edge_vals[:, 1]) / (c_edge_vals[:, 0] - c_edge_vals[:, 1])
        w1 = (-c + c_edge_vals[:, 0]) / (c_edge_vals[:, 0] - c_edge_vals[:, 1])

        # Calculate points linearly interpolated from vertices
        points = mesh.vertices[mesh.edges_unique[edge_inds]]
        points = points[:, 0] * w0[:, None] + points[:, 1] * w1[:, None]

        # Determine adjacency
        c_edges_in_faces = edge_inds[mesh.faces_unique_edges]
        c_faces = np.any(c_edges_in_faces, axis=1)
        edges_in_c_faces = mesh.faces_unique_edges[c_faces]

        if len(edges_in_c_faces) == 0:
            print("No contours at f=", c)
            continue

        # Each element of this array corresponds to a face containing the contour
        # The two values in the element are the edges (indices) adjacent to the face
        c_edges_in_c_faces = np.array(
            [a[b] for a, b in zip(edges_in_c_faces, c_edges_in_faces[c_faces])]
        )

        # Contour edges (indices pointing to the unique edge list)
        c_edge_inds = list(np.flatnonzero(edge_inds))

        # Indices of faces used in contour
        c_face_inds = np.flatnonzero(c_faces)

        # Check gradient of stream function in first (could be any) triangle
        c_face_gradient = g[:, c_face_inds[0]]

        # Vector between two edges in first face used for contour
        vec = (
            points[c_edge_inds.index(c_edges_in_c_faces[0, 0])]
            - points[c_edge_inds.index(c_edges_in_c_faces[0, 1])]
        )

        # Create loop such that it is in the same direction as gradient
        if c_face_gradient.dot(vec) >= 0:
            # Init loop variables
            key = c_edges_in_c_faces[0, 1]
            val = c_edges_in_c_faces[0, 0]
        else:
            # Init loop variables
            key = c_edges_in_c_faces[0, 0]
            val = c_edges_in_c_faces[0, 1]

        ii = 0
        sorted_inds = []
        kmax = len(c_edge_inds)
        k = 0

        # Loop over c_edges by essentially solving a linked list from c_edges_in_c_faces
        while k < kmax:
            sorted_inds.append(c_edge_inds.index(val))
            c_edges_in_c_faces[ii] = -1
            ii, jj = np.nonzero(c_edges_in_c_faces == key)

            if len(ii) == 0:
                # Next edge not found in the adjacency list, contour must be closed now
                # OR the loop is not closed, test this:
                ii_temp, jj_temp = np.nonzero(c_edges_in_c_faces == val)
                if len(ii_temp) > 0:
                    # Proceed to another direction (direction does not matter)
                    ii = ii_temp
                    jj = jj_temp
                else:
                    # Sort points containing contours by adjacency of the edges
                    # and append to contour_polys

                    contour_polys.append(points[sorted_inds])
                    contour_values.append(c)

                    # Break the loop if all edges have been visited
                    if np.all(c_edges_in_c_faces == -1):
                        break

                    # Else find a starting point in another contour at the same level
                    sorted_inds = []

                    ii = np.flatnonzero(c_edges_in_c_faces[:, 0] >= 0)[0]
                    jj = 0

                    # Compute gradient in face
                    c_face_gradient = g[:, c_face_inds[ii]]

                    # once again, check winding direction
                    vec = (
                        points[c_edge_inds.index(c_edges_in_c_faces[ii, 0])]
                        - points[c_edge_inds.index(c_edges_in_c_faces[ii, 1])]
                    )

                    if c_face_gradient.dot(vec) >= 0:
                        jj = 0
                    else:
                        jj = 1

            else:
                # Edge found
                ii = ii[0]
                jj = jj[0]

            # Update key and value
            val = c_edges_in_c_faces[ii, jj]
            key = c_edges_in_c_faces[ii, (jj + 1) % 2]

            k += 1

            if k == kmax:
                raise RuntimeWarning(
                    "Something wrong with the contours, number of max iterations exceeded"
                )
    if return_values:
        return contour_polys, contour_values

    return contour_polys


def simplify_contour(c, min_edge=1e-3, angle_threshold=2e-2, smooth=True):
    """
    Simplifies contours by merging small (short) segments and
    with only a small angle difference.

    Optionally applies smoothing to contour shape.

    Parameters
    ----------
    c: list
        List of polygons describing closed loops.
    min_edge: float
        Minimum edge length. Edges shorter than this are merged.
    angle_threshold: float
        Minimum angle. Edges with smaller angle differences are merged.
    smooth: bool
        If True, apply smoothing to the polygon shapes.

    Returns
    -------
    c: list
        Modified list of polygons

    """
    # Remove small edges by threshold
    vals = [np.ones(c.shape[0]), -np.ones(c.shape[0]), np.ones(c.shape[0])]
    D = spdiags(vals, [1, 0, -c.shape[0] + 1], c.shape[0], c.shape[0])
    edges = D @ c
    c = c[np.linalg.norm(edges, axis=1) > min_edge]
    if len(c) == 0:
        return None

    # Remove nodes on straight lines
    D = spdiags(vals, [1, 0, -c.shape[0] + 1], c.shape[0], c.shape[0])
    H = spdiags(1 / np.linalg.norm(D @ c, axis=1), 0, c.shape[0], c.shape[0])
    DD = H @ D
    c = c[np.linalg.norm(D.T @ DD @ c, axis=-1) > angle_threshold]

    if smooth:
        D = spdiags(vals, [1, 0, -c.shape[0] + 1], c.shape[0], c.shape[0])
        H = spdiags(1 / np.linalg.norm(D @ c, axis=1), 0, c.shape[0], c.shape[0])
        DD = H @ D
        lengths = np.linalg.norm(D @ c, axis=1)
        lengths = 0.5 * abs(D.T) @ lengths  # Mean of edges
        #            c = c - 0.2*lengths[:,None]*(D.T @ DD @ c)
        Nc = c.shape[0]
        c = spsolve(speye(Nc, Nc) + 1.0 * spdiags(lengths, 0, Nc, Nc) @ (D.T @ DD), c)

    return c
