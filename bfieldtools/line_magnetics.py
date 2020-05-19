"""
Functions for working with current polylines, e.g. for calculating the magnetic field and potentials as well as the inductance

"""

__all__ = [
    "magnetic_field",
    "mutual_inductance",
    "scalar_potential",
    "self_inductance",
    "vector_potential",
]

import numpy as np
from .integrals import omega
import trimesh


def cross(r1, r2):
    """ Cross product without overhead
    """
    result = np.zeros(r1.shape)
    result[0] = r1[1] * r2[2] - r1[2] * r2[1]
    result[1] = r1[2] * r2[0] - r1[0] * r2[2]
    result[2] = r1[0] * r2[1] - r1[1] * r2[0]
    return result


def magnetic_field(vertices, points):
    """ Compute B-field of a segmented line current.
        See: Compact expressions for the Biot–Savart fields of a filamentary segments
        by Hanson & Hirshman: https://doi.org/10.1063/1.1507589


        Parameters
        ----------

        vertices: (N_line, 3) array
            Vertices of the line with N_line-1 segments
            The first and last vertices should be the same to close the loop.
        points:   (N_points, 3) array
            Magnetic field evaluation points

        Returns
        -------
        bfield: (N_points, 3) array
            Magnetic field at evaluation points

    """
    field = np.zeros(points.T.shape)
    for i in range(len(vertices) - 1):
        r1 = vertices[i]
        r2 = vertices[i + 1]

        # Vectors between vertices and field points
        a1 = points.T - r1.reshape(3, 1)
        a2 = points.T - r2.reshape(3, 1)

        # Direction of the field
        f = cross(a1, a2)

        # Vector lengths
        d1 = np.sqrt(np.sum(a1 ** 2, axis=0))
        d2 = np.sqrt(np.sum(a2 ** 2, axis=0))

        # Normalize direction field and divide by cylindrical distance
        f *= (d1 + d2) / (d1 * d2 * (d1 * d2 + np.sum(a1 * a2, axis=0)))

        field = field + f

    return field.T * 1e-7


def vector_potential(vertices, points, reg=1e-12, symmetrize=True):
    """ Compute vector potential of a segmented line currents.
        Based on straightforward integration of 1/r potential over a line
        i.e. the gamma0 integral

        See: Compact expressions for the Biot–Savart fields of a filamentary segments
        by Hanson & Hirshman: https://doi.org/10.1063/1.1507589


        Parameters
        ----------
        vertices: (N_line, 3) array
            Vertices of the line with N_line-1 segments
            The first and last vertices should be the same to close the loop.
        points: (N_points, 3) array
            Evaluation points

        Returns
        -------
        A: array (Npoints, 3)
            Vector potential

    """

    segments = vertices[1:] - vertices[:-1]
    RR = vertices[:, None, :] - points[None, :, :]
    dotprods2 = np.sum(RR[1:] * segments[..., None, :], axis=-1)
    dotprods1 = np.sum(RR[:-1] * segments[..., None, :], axis=-1)
    ss = np.linalg.norm(segments, axis=-1)
    segments /= ss[..., None]
    rr = np.linalg.norm(RR, axis=-1)

    # Regularize s.t. neither the denominator or the numerator can be zero
    # Avoid numerical issues directly at the edge
    res = np.log(
        (rr[1:] * ss[..., None] + dotprods2 + reg)
        / (rr[:-1] * ss[..., None] + dotprods1 + reg)
    )

    # Symmetrize the result since on the negative extension of the edge
    # there's division of two small values resulting numerical instabilities
    # (also incompatible with adding the reg value)
    if symmetrize:
        res2 = -np.log(
            (rr[1:] * ss[..., None] - dotprods2 + reg)
            / (rr[:-1] * ss[..., None] - dotprods1 + reg)
        )
        res = np.where(dotprods1 + dotprods2 > 0, res, res2)

    return 1e-7 * np.sum(res[..., None] * segments[..., None, :], axis=0)


def scalar_potential(vertices, points):
    """ Computes the scalar magnetic potential of a segmented current loop at given points.
        This is equal to the solid angle spanned by the loop (polygon), times a constant.
        The first and last vertices should be the same to close the loop.


        Parameters
        ----------
        vertices: (N_line, 3) array
            Vertices of the line with N_line-1 segments
            The first and last vertices should be the same to close the loop.
        points: (N_points, 3) array
            Evaluation points

        Returns
        -------
        U: array (Npoints, )
            Scalar magnetic potential

    """

    N_verts = len(vertices)

    # VERTEX MASS CENTRE
    mass_center = np.mean(vertices, axis=0)

    vertices = np.vstack((vertices, mass_center))

    # CREATE TRIANGLE FAN
    faces = np.full(shape=(N_verts - 1, 3), fill_value=np.nan, dtype=int)

    for i in range(N_verts - 1):
        faces[i] = np.array([i, i + 1, N_verts])

    R1 = vertices[faces]
    R2 = points

    RR = R2[:, None, None, :] - R1[None, :, :, :]

    # COMPUTE SOLID ANGLE, DIVIDE BY 4*PI*mu0
    return np.sum(omega(RR), axis=1) * 1e-7


def mutual_inductance(path1, path2, Nquad=2, radius=1e-3):
    """ Compute magnetic flux created by a segmented line current loops
        (path1) on a another closed loop of segmented current
        (path2). The other loop(s) is numerically integrated.

        In other words, calculate mutual inductance of the current loops.

        **SELF INDUCTANCE OF A LOOP IS BASED ON ROUND-WIRE APPROXIMATION**
        
        Parameters
        ----------
        path1: trimesh.Path3D-object
        
        path2: trimesh.Path3D-object

        Nquad: int, the number of quadrature points on line segments
        
        radius: float
            radius of the wire, only used for self-inductance calcultion
            self-inductance is dependent on the wire cross section
        Returns
        -------
         flux in the other loop generated by loops (Nloops,)


    """
    # Calculate quadrature points linearly spaced on the line segments
    # of the other loop
    t = np.linspace(0, 1, Nquad + 1)

    fluxes = np.zeros((len(path1.entities), len(path2.entities)))
    for j, loop2 in enumerate(path2.entities):
        vertices_other = path2.vertices[loop2.nodes, :]
        # sides = vertices_other[1:] - vertices_other[:-1]
        sides = vertices_other[:, 1] - vertices_other[:, 0]
        segments = (t[1:, None, None] - t[:-1, None, None]) * sides

        tq = 0.5 * (t[1:] + t[:-1])  # Nquad points on the segments

        # Nquad, Nsegments, 3
        points = vertices_other[:, 0] + tq[:, None, None] * sides
        shape = points.shape

        for i, loop1 in enumerate(path1.entities):
            if (path1 is path2) and i == j:
                fluxes[i, j] = self_inductance(path1, loop1, radius)
            else:
                vertices = path1.vertices[loop1.points]
                a = vector_potential(vertices, points.reshape(-1, 3))

                # Nquad, Nsegments, 3
                a = a.reshape(shape)

                # Take dot product between vector potential and the line segments
                # corresponding to each quadrature points (axis=0) and sum over the
                # segments (axis=1) and quadrature points on each segment (axis=2)

                fluxes[i, j] = np.sum(a * segments, axis=(0, 1, 2))

    return fluxes


def self_inductance(path, loop, radius=1e-3):
    """
    Calculate self inductance based on round-wire approximation
    http://www.thompsonrd.com/induct2.pdf
    section 5.
    

    Parameters
    ----------
    path : path3D object
    loop : entity of the path3D object
        entity describing the single loop
    radius : the radius of the wire
        DESCRIPTION. The default is 1e-3.

    Returns
    -------
    self_ind : float

    """
    mu0_over2pi = (1e-7) * 2
    # Perimeter
    p = loop.length(path.vertices)

    # CREATE TRIANGLE FAN for area calculation
    vertices = path.vertices[loop.points]
    N_verts = len(vertices)
    mass_center = np.mean(vertices, axis=0)
    vertices = np.vstack((vertices, mass_center))
    faces = np.full(shape=(N_verts - 1, 3), fill_value=np.nan, dtype=int)
    for i in range(N_verts - 1):
        faces[i] = np.array([i, i + 1, N_verts])
    mesh = trimesh.Trimesh(vertices, faces)
    # Area
    A = mesh.area

    # This is an approximation for round wire
    self_ind = mu0_over2pi * p * (np.log(2 * A / (p * radius)) + 0.25)

    return self_ind
