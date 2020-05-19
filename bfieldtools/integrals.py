"""

Analytic integral for vectorized field / potential computation

"""

__all__ = [
    "c_coeffs",
    "d_distance",
    "gamma0",
    "omega",
    "potential_dipoles",
    "potential_vertex_dipoles",
    "triangle_potential_approx",
    "triangle_potential_dipole_linear",
    "triangle_potential_uniform",
    "x_distance",
    "x_distance2",
]

import numpy as np


def determinant(a):
    """ Faster determinant for the two last dimensions of 'a'
    """
    det = a[..., 0, 0] * (a[..., 1, 1] * a[..., 2, 2] - a[..., 2, 1] * a[..., 1, 2])
    det += a[..., 0, 1] * (a[..., 1, 2] * a[..., 2, 0] - a[..., 2, 2] * a[..., 1, 0])
    det += a[..., 0, 2] * (a[..., 1, 0] * a[..., 2, 1] - a[..., 2, 0] * a[..., 1, 1])
    return det


def norm(vecs):
    """ Faster vector norm for the last dimension of 'vecs'
    """
    return np.sqrt(np.einsum("...i,...i", vecs, vecs))


def cross(r1, r2):
    """ Cross product without overhead for the last dimensions of 'r1' and 'r2'
    """
    result = np.zeros(r1.shape)
    result[..., 0] = r1[..., 1] * r2[..., 2] - r1[..., 2] * r2[..., 1]
    result[..., 1] = r1[..., 2] * r2[..., 0] - r1[..., 0] * r2[..., 2]
    result[..., 2] = r1[..., 0] * r2[..., 1] - r1[..., 1] * r2[..., 0]
    return result


def gamma0(R, reg=1e-13, symmetrize=True):
    """ 1/r integrals over the edges of a triangle called gamma_0
        (line charge potentials).

        **NOTE: MAY NOT BE VERY PRECISE FOR POINTS DIRECTLY AT TRIANGLE
        EDGES.**

        Parameters
        ----------
        R : ndarray (..., N_triverts, xyz)
            displacement vectors (r-r') between Neval evaluation points (r)
            and the 3 vertices of the Ntri triangles/triangle.
        reg: float, a small value added to the arguments of the logarithm,
             regularizes the values very close to the line segments
        symmetrize: recalculates the result for by mirroring
                    the evaluation points with respect the line segment
                    mid point to get rid off the badly behaving points on the
                    negative extension of the line segment


        Returns
        -------
        res: array (Neval, Nverts)
            The analytic integrals for each vertex/edge

    """
    edges = np.roll(R[0], 2, -2) - np.roll(R[0], 1, -2)
    dotprods1 = np.einsum("...i,...i", np.roll(R, 1, -2), edges)
    dotprods2 = np.einsum("...i,...i", np.roll(R, 2, -2), edges)
    en = norm(edges)
    del edges
    n = norm(R)
    # Regularize s.t. neither the denominator or the numerator can be zero
    # Avoid numerical issues directly at the edge
    nn1 = np.roll(n, 1, -1) * en
    nn2 = np.roll(n, 2, -1) * en
    res = np.log((nn1 + dotprods1 + reg) / (nn2 + dotprods2 + reg))

    # Symmetrize the result since on the negative extension of the edge
    # there's division of two small values resulting numerical instabilities
    # (also incompatible with adding the reg value)
    if symmetrize:
        mask = ((np.abs(dotprods1 + nn1)) < 1e-12) * (dotprods1 + dotprods2 < 0)
        res[mask] = -np.log(
            (nn1[mask] - dotprods1[mask] + reg) / (nn2[mask] - dotprods2[mask] + reg)
        )

    res /= en
    return -res


def omega(R):
    """ Calculate the solid angle of a triangles

        see
        A. Van Oosterom and J. Strackee
        IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING,
        VOL. BME-30, NO. 2, 1983

        Parameters
        ----------
        R : ndarray (Neval, (Ntri), N_triverts, xyz)
            displacement vectors (r-r') of Ntri triangles
            and Neval evaluation points for the 3 vertices
            of the triangles/triangle.

            The shape of R can any with the constraint that
            the last dimenion corrsponds to coordinates (x, y, z) and the
            second last dimension to triangle vertices (vert1, vert2, vert3)

        Returns
        -------
        sa: (Neval, (Ntri))
            Solid angles of subtened by triangles at evaluation points
    """
    # Distances
    d = norm(R)
    # Scalar triple products
    stp = determinant(R)
    # Denominator
    denom = np.prod(d, axis=-1)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        #        denom += np.sum(R[..., i, :]*R[..., j, :], axis=-1)*d[..., k]
        denom += np.einsum("...i,...i,...", R[..., i, :], R[..., j, :], d[..., k])
    # Solid angles
    sa = -2 * np.arctan2(stp, denom)
    return sa


def x_distance(R, tn, ta=None):
    """ Signed distances in the triangle planes from the opposite
        edge towards the node for all evaluation points in R

        The distances are normalized to one at the node if areas are given
        The distances are multiplied by the edge lenght if areass are None

        Parameters:

            R: ndarray (... Ntri, Nverts, xyz)
                displacement vectors (coordinates)
            tn: ndarray (Ntri, 3)
                triangle normals
            ta: ndarray (Ntri)
                triangle areas
                if None, normalizization with double area is not carried out

        returns:
            ndaarray (..., Ntri, N_triverts (3)), distance in the triangle plane
    """
    edges = np.roll(R[0], 2, -2) - np.roll(R[0], 1, -2)
    if ta is not None:
        edges /= 2 * ta[:, None, None]
    edges = -cross(edges, tn[:, None, :])
    return np.einsum("...k,...k->...", np.roll(R, 1, -2), edges)


def x_distance2(mesh):
    """ Signed distances in the triangle planes from the opposite
        edge towards the node for all evalution points in R
    """
    # TODO: with gradient, needs mesh info
    pass


def d_distance(R, tn):
    """ Signed distance from the triangle plane for each triangle

        Parameters:

            R: ndarray (... Ntri, Nverts, xyz)
                displacement vectors (coordinates)
            tn: ndarray (Ntri, 3)
                triangle normals

        Returns:

            ndarray (..., Ntri, N_triverts (3)) of signed distances
    """
    return np.einsum("...ki,ki->...k", np.take(R, 0, -2), tn)


def c_coeffs(R, ta):
    """ Cotan-coeffs

        Parameters:

            R: ndarray (... Ntri, Nverts, xyz)
                displacement vectors (coordinates)
            ta: ndarray (Ntri)
                triangle areas


        Returns:

            ndarray (..., Ntri, N_triverts (3))
    """
    edges = np.roll(R[0], 2, -2) - np.roll(R[0], 1, -2)
    return np.einsum("...ik,...jk->...ij", edges, edges / (2 * ta[:, None, None]))


def triangle_potential_uniform(R, tn, planar=False):
    """ 1/r potential of a uniform triangle

        for original derivation see
        A. S. Ferguson, Xu Zhang and G. Stroink,
        "A complete linear discretization for calculating the magnetic field
        using the boundary element method,"
        in IEEE Transactions on Biomedical Engineering,
        vol. 41, no. 5, pp. 455-460, May 1994.
        doi: 10.1109/10.293220


        Parameters
        ----------

        R : (Neval, (Ntri), 3, 3) array
            Displacement vectors (Neval, (Ntri), Ntri_verts, xyz)
        tn : ((Ntri), 3) array
            Triangle normals (Ntri, dir)
        planar: boolean
            If True, assume all the triangles and the evaluation points
                    are on the same plane (for speed), leaves out the
                    omega term

        Returns
        -------
        result: result:  ndarray (Neval, (Ntri))
            Resultant 1/r potential for each triangle (Ntri)
            at the field evaluation points (Neval)

    """
    x = x_distance(R, tn, None)
    result = np.einsum("...i,...i", gamma0(R), x)
    if not planar:
        result += d_distance(R, tn) * omega(R)
    return result


def triangle_potential_approx(Rcenters, ta, reg=1e-12):
    """ 1/r potential of a uniform triangle using centroid approximation

        Calculates 1/R potentials for triangle centroids
        (The singularity at the centroid is handled with the very small
        reg value, but anyway the values close to the centroid are inexact)

        Parameters
        ----------
        Rcenters : (N, (Ntri), 3) array
            Displacement vectors (Neval, Ntri, xyz)
            from triangle centers
        ta : (Ntri) array
            Triangle areas

        reg: float
            Regularization value used in approximation

        Returns
        -------
        result: result:  ndarray (...., Ntri, Ntri_verts)
            Resultant 1/r potential for each node (Ntri_verts)
            in each triangle (Ntri) in the displacement vectors R

    """
    result = ta / (norm(Rcenters) + reg)
    return result


def potential_dipoles(R, face_normals, face_areas):
    """ Approximate the potential of linearly varying dipole density by
        by dipoles at each face

    Parameters
            R : ndarray (Neval, Ntri, Ntri_verts, N_xyz)
                Displacement vectors
            face_normals: ndarray (Ntri, 3)
                normals for each triangle
            face_areas: ndarray (Ntri,)
                areas for each triangle

    Return
        Potential approximation for vertex in each face
        pot: (Neval, Ntri, Ntriverts)
    """
    nn = face_normals
    # Calculate quadrature points corresponding to linear shape functions (Ok?)
    weights = np.array([[0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5]])
    #    weights = np.eye(3)
    #    weights = np.ones((3,3))/3
    # Combine vertices for quadrature points
    Rquad = np.einsum("...ij,ik->...kj", R, weights)
    pot = np.einsum("ik, ...ijk->...ij", nn, Rquad) / (norm(Rquad) ** 3)
    pot = pot * (face_areas[:, None] / 3)

    return pot


def potential_vertex_dipoles(R, vertex_normals, vertex_areas):
    """ Approximate the potential of linearly varying dipole density by
        by dipoles at each vertex

    Parameters
            R : ndarray (Neval, Nvertex, N_xyz)
                Displacement vectors
            vertex_normals: ndarray (Nvertex, 3)
                normals for each triangle
            vertex_areas: ndarray (Nvertex,)
                areas for each triangle

    Return
        Potential approximation for vertex in each face
        pot: (Neval, Ntri, Ntriverts)
    """
    nn = vertex_normals
    pot = np.einsum("ik, lik->li", nn, R) / (norm(R) ** 3)
    pot *= vertex_areas

    return pot


def triangle_potential_dipole_linear(R, tn, ta):
    """ Potential of dipolar density with magnitude of a
        linear shape function on a triangle, "omega_i" in de Munck's paper

        for the original derivation, see:
        J. C. de Munck, "A linear discretization of the volume mesh_conductor
        boundary integral equation using analytically integrated elements
        (electrophysiology application),"
        in IEEE Transactions on Biomedical Engineering,
        vol. 39, no. 9, pp. 986-990, Sept. 1992.
        doi: 10.1109/10.256433

        Parameters
        ----------

        R : (..., Ntri, 3, 3) array
            Displacement vectors (...., Ntri, Ntri_verts, xyz)
        tn : ((Ntri), 3) array
            Triangle normals (Ntri, dir)
        ta : (Ntri), array
            Triangle areas (Ntri, dir)

        Returns
        -------
        result:  ndarray (...., Ntri, Ntri_verts)
            Resultant dipolar potential for each shape functions (Ntri_verts)
            in each triangle (Ntri) at the points
            corresponding to displacement vectors in R

    """

    result = np.einsum(
        "...i,...ij,...->...j",
        gamma0(R),
        c_coeffs(R, ta),
        d_distance(R, tn),
        optimize=True,
    )
    x_dists = x_distance(R, tn, ta)
    result -= x_dists * omega(R)[..., :, None]

    return result
