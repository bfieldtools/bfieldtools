"""
Analytic integral for vectorized field / potential computation

"""

import numpy as np


def determinant(a):
    det = a[..., 0, 0] * (a[..., 1, 1] * a[..., 2, 2] - a[..., 2, 1] * a[..., 1, 2])
    det += a[..., 0, 1] * (a[..., 1, 2] * a[..., 2, 0] - a[..., 2, 2] * a[..., 1, 0])
    det += a[..., 0, 2] * (a[..., 1, 0] * a[..., 2, 1] - a[..., 2, 0] * a[..., 1, 1])
    return det


def norm(vecs):
    return np.sqrt(np.einsum("...i,...i", vecs, vecs))


def cross(r1, r2):
    """ Cross product without overhead
    """
    result = np.zeros(r1.shape)
    result[..., 0] = r1[..., 1] * r2[..., 2] - r1[..., 2] * r2[..., 1]
    result[..., 1] = r1[..., 2] * r2[..., 0] - r1[..., 0] * r2[..., 2]
    result[..., 2] = r1[..., 0] * r2[..., 1] - r1[..., 1] * r2[..., 0]
    return result


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
    edges = np.roll(R[0], 1, -2) - np.roll(R[0], 2, -2)
    #    dotprods1 = np.sum(np.roll(R, 1, -2)*edges, axis=-1)
    #    dotprods2 = np.sum(np.roll(R, 2, -2)*edges, axis=-1)
    dotprods1 = np.einsum("...i,...i", np.roll(R, 1, -2), edges)
    dotprods2 = np.einsum("...i,...i", np.roll(R, 2, -2), edges)
    en = norm(edges)
    del edges
    n = norm(R)
    # Regularize s.t. neither the denominator or the numerator can be zero
    # Avoid numerical issues directly at the edge
    nn1 = np.roll(n, 2, -1) * en
    nn2 = np.roll(n, 1, -1) * en
    res = np.log((nn1 + dotprods2 + reg) / (nn2 + dotprods1 + reg))

    # Symmetrize the result since on the negative extension of the edge
    # there's division of two small values resulting numerical instabilities
    # (also incompatible with adding the reg value)
    if symmetrize:
        res2 = -np.log((nn1 - dotprods2 + reg) / (nn2 - dotprods1 + reg))
        res = np.where(dotprods1 + dotprods2 > 0, res, res2)
    res /= en
    return -res  # TODO: there should be minus, since we want this to be positive


def omega(R):
    """ Calculate the solid angle of a triangles

        see
        A. Van Oosterom and J. Strackee
        IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING,
        VOL. BME-30, NO. 2, 1983

        Parameters
        ----------
        R : array of points (Neval, (Ntri), Nverts, xyz)
            Points correspond to relative coordinates (x,y,z) of
            N triangles/evaluation points for
            the 3 corners of the triangles/triangle.

            Neval can be number of evaluation points for the same triangle
            or number of triangles for the same evaluation points

        Returns
        -------
        sa: (Neval, (Ntri))
            Solid angles of triangle(s) at evaluation points
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
    sa = 2 * np.arctan2(stp, denom)
    return sa


def triangle_potential_uniform(R, tn, planar=False):
    """ 1/r potential of a uniform triangle

        see
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
            If True, use planar geometry assumption for speed

        Returns
        -------
        result: result:  ndarray (Neval, (Ntri))
            Resultant 1/r potential for each triangle (Ntri)
            at the field evaluation points (Neval)

    """
    if len(R.shape) > 3:
        tn_ax = tn[:, None, :]
    else:
        tn_ax = tn
    summands = np.sum(
        tn_ax * np.cross(np.roll(R, 1, -2), np.roll(R, 2, -2), axis=-1), axis=-1
    )
    #    summands = -gamma0(R)*np.sum(tn_ax*np.cross(np.roll(R, 1, -2),
    #                                               np.roll(R, 2, -2), axis=-1), axis=-1)
    result = np.einsum("...i,...i", -gamma0(R), summands)
    if not planar:
        #        csigned = np.sum(np.take(R, 0, -2)*tn, axis=-1)
        csigned = np.einsum("...i,...i", np.take(R, 0, -2), tn)
        result -= csigned * omega(R)
    else:
        print("Assuming all the triangles are in the same plane!")
    #        result = np.sum(summands, axis=-1)
    return result


def triangle_potential_approx(R, ta, reg=1e-12):
    """ 1/r potential of a uniform triangle using centroid approximation

        Calculates 1/R potentials for triangle centroids
        (The singularity at the centroid is handled with the very small
        reg value, but anyway the values close to the centroid are inexact)

        Parameters
        ----------
        R : (N, (Ntri), 3, 3) array
            Displacement vectors (Neval, ...., Ntri_verts, xyz)
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
    result = 1 / (norm(np.mean(R, axis=-2)) + reg) * ta
    return result


def potential_dipoles(R, face_normals, face_areas):
    """ Approximate the potential of linearly varying dipole density by
        by dipoles at each face

    Parameters
            R : (Neval, Ntri, 3, 3) array
                Displacement vectors (Neval, Ntri, Ntri_verts, xyz)
            face_normals: (Ntri, 3) normals of each fame
            face_areas

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


def triangle_potential_dipole_linear(R, tn, ta, planar=False):
    """ Potential of dipolar density with magnitude of a
        linear shape function on a triangle, "omega_i" in de Munck's paper

        see
        J. C. de Munck, "A linear discretization of the volume conductor
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
        planar: boolean
            If True, use planar geometry assumption for speed

        Returns
        -------
        result:  ndarray (...., Ntri, Ntri_verts)
            Resultant dipolar potential for each shape functions (Ntri_verts)
            in each triangle (Ntri) at the points
            corresponding to displacement vectors in R

    """
    if len(R.shape) > 3:
        tn_ax = tn[:, None, :]
    else:
        tn_ax = tn
    # Volumes of tetrahedron between field evaluation point and the triangle
    #    det = np.sum(np.cross(np.roll(R, 2, -2),
    #                          np.roll(R, 1, -2), axis=-1)*R, axis=-1)
    det = determinant(R)
    # Edges opposite to the nodes
    edges = np.roll(R[0], 1, -2) - np.roll(R[0], 2, -2)
    # Latter part of omega_i integral in de Munck
    #    result = np.sum(np.sum(gamma0(R)[..., None]*edges, axis=-2)[...,None,:]*edges, axis=-1)
    result = np.einsum(
        "...i,...ik,...jk,...->...j",
        gamma0(R),
        edges,
        edges,
        det / (2 * ta),
        optimize=True,
    )
    #    result *= (det/(2*ta))[..., None] # TODO: IS DET SIGN OK?
    if not planar:
        # First part of the integral
        # Note: tn normalized version of n-vector in de Munck
        lin_coeffs = np.sum(
            tn_ax * cross(np.roll(R, 2, -2), np.roll(R, 1, -2)), axis=-1
        )
        result += lin_coeffs * omega(R)[..., :, None]
    else:
        print("Assuming all the triangles are in the same plane!")
    return result / (2 * ta[:, None])
