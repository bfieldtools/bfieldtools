'''
Functions for working with current line segments

'''


import numpy as np
from .integrals import omega


def cross(r1, r2):
    """ Cross product without overhead
    """
    result = np.zeros(r1.shape)
    result[0] = r1[1] * r2[2] - r1[2] * r2[1]
    result[1] = r1[2] * r2[0] - r1[0] * r2[2]
    result[2] = r1[0] * r2[1] - r1[1] * r2[0]
    return result


def magnetic_field2(vertices, points):
    """ Compute b field of a segmented line current.
        See:
        Compact expressions for the Biot–Savart fields of a filamentary segments
        by Hanson & Hirshman


        Parameters
        ----------

        vertices: (N_line, 3) array
            Vertices of the line with N_line-1 segments
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
        d1 = np.sqrt(np.sum(a1**2, axis=0))
        d2 = np.sqrt(np.sum(a2**2, axis=0))

        # Normalize direction field and divide by cylindrical distance
        f *= (d1 + d2)/(d1*d2*(d1*d2 + np.sum(a1*a2, axis=0)))

        field = field + f

    return field.T * 1e-7


def magnetic_field(vertices, points):
    """ Compute b field of a segmented line current.
        This calculation is based on integration by Griffiths
        on page 217 (3rd edition)

        Parameters
        ----------

        vertices: (N_line, 3) array
            Vertices of the line with N_line-1 segments
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
        d = ((r2 - r1)/((r1 - r2)**2).sum()).reshape(3, 1)

        # Vectors between vertices and field points
        a1 = points.T - r1.reshape(3, 1)
        a2 = points.T - r2.reshape(3, 1)

        # Direction of the field
        f = cross(a1, d)

        # Sine factor
        sinefactor = (d * a2).sum(axis=0) / np.sqrt((a2**2).sum(axis=0))
        sinefactor = sinefactor - (d * a1).sum(axis=0) / np.sqrt((a1**2).sum(axis=0))

        # Normalize direction field and divide by cylindrical distance
        s2 = (f**2).sum(axis=0)
        s2[s2 == 0] = 1e-12 # Regularize for points directly at the
                            # continuation of the line segment
        f *= (sinefactor / s2)

        field = field + f

    return field.T * 1e-7


def vector_potential(vertices, points,
                     reg=1e-12, symmetrize=True):
    """ Compute vector potential of a segmented line currents.
        Based on straightforward integration of 1/r potential over a line
        i.e. the gamma0 integral

        See:
            Compact expressions for the Biot–Savart fields of a filamentary segments
            by Hanson & Hirshman


        Parameters
        ----------
        vertices: (N_line, 3) array
            Vertices of the line with N_line-1 segments
        points: (N_points, 3) array
            Magnetic field evaluation points
        Returns
        -------
            Vector potential (Nloops, Npoints, 3)

    """
 
    loops = np.array([np.arange(len(vertices))])

    loops2 = np.roll(loops, -1, 1)
    loops1 = loops
    segments = vertices[loops2] - vertices[loops1]
    RR = vertices[:, None, :] - points[None, :, :]
    dotprods2 = np.sum(RR[loops2] * segments[..., None, :], axis=-1)
    dotprods1 = np.sum(RR[loops1] * segments[..., None, :], axis=-1)
    ss = np.linalg.norm(segments, axis=-1)
    segments /= ss[..., None]
    rr = np.linalg.norm(RR, axis=-1)

    # Regularize s.t. neither the denominator or the numerator can be zero
    # Avoid numerical issues directly at the edge
    res = np.log((rr[loops2] * ss[..., None] + dotprods2 + reg)
                 / (rr[loops1] * ss[..., None] + dotprods1 + reg))

    # Symmetrize the result since on the negative extension of the edge
    # there's division of two small values resulting numerical instabilities
    # (also incompatible with adding the reg value)
    if symmetrize:
        res2 = -np.log((rr[loops2] * ss[..., None] - dotprods2 + reg)
                       / (rr[loops1] * ss[..., None] - dotprods1 + reg))
        res = np.where(dotprods1 + dotprods2 > 0, res, res2)

    return 1e-7 * np.sum(res[..., None] * segments[..., None, :], axis=1)


def scalar_potential(vertices, points):
    '''
    Computes the scalar magnetic potential of a segmented current loop at given points.
    This is equal to the solid angle spanned by the loop (polygon), times a constant.
    The first and last vertices are connected to close the loop.
    Parameters
    ----------
    vertices: (N_line, 3) array
        Vertices of the line with N_line-1 segments
    points: (N_points, 3) array
        Magnetic field evaluation points

    Returns
    -------
    Scalar magnetic potential (Npoints, )

    '''

    N_verts = len(vertices)

    #VERTEX MASS CENTRE
    mass_center = np.mean(vertices, axis=0)

    vertices = np.vstack((vertices, mass_center))

    #CREATE TRIANGLE FAN
    faces = np.full(shape=(N_verts, 3), fill_value=np.nan, dtype=int)

    for i in range(N_verts):
        faces[i] = np.array([i, (i+1)%N_verts, N_verts])

    R1 = vertices[faces]
    R2 = points

    RR = R2[:, None, None, :] - R1[None, :, :, :]

    #COMPUTE SOLID ANGLE, DIVIDE BY 4*PI*mu0
    return np.sum(omega(RR), axis=1) * 1e-7


def flux(vertices, loops, vertices_other, Nquad=2):
    """ Compute magnetic flux created by a segmented line current loops
        (vertices, loops) on a another closed loop of segmented current
        (vertices_other). The other loop is numerically integrated.

        In other words, calculate mutual inductance of the current loops.

        NOT SUITABLE for calculating the self-flux, i.e., self inductance

        Parameters
        ----------
        vertices:
            all vertices in segmented loops generating the flux

        loops:
            list of indices defining closed loops of vertices, if
            None use all vertices. All loops must have the same
            number of indices (this could be changed in future)
            Example: Giving array of 4 vertices, the loops can be
            defined as loops = np.array([[0,1,2,3]])

        vertices_other:
            vertices in the loop receiving the flux


        Returns
        -------
         flux in the other loop generated by loops (Nloops,)


    """
    # Calculate quadrature points linearly spaced on the line segments
    # of the other loop
    t = np.linspace(0, 1, Nquad + 1)
    sides = vertices_other[1:] - vertices_other[:-1]
    segments = (t[1:, None, None] - t[:-1, None, None]) * sides

    t = 0.5 * (t[1:] + t[:-1]) # Nquad points on the segments

    # Nquad, Nsegments, 3
    points = vertices_other[:-1] + t[:, None, None] * sides
    shape = points.shape

    # Nloops, Nquad*Nsegments, 3
    a = vector_potential(vertices, points.reshape(-1, 3), loops)

    # Nloops, Nquad, Nsegments, 3
    a = a.reshape(a.shape[0:1] + shape)

    # Take dot product between vector potential and the line segments
    # corresponding to each quadrature points (axis=3) and sum over the
    # segements (axis=2) and quadrature points on each segment (axis=1)

    return np.sum(a * segments, axis=(1, 2, 3))


#
#if __name__ == "__main__":
#    """
#    Plot field of a circular current path
#
#    """
#    x = np.linspace(-1, 1, 200)
#    Ntheta = 5
#    theta = np.linspace(0, 2 * np.pi, Ntheta)
#    vertices = np.zeros((Ntheta, 3), dtype=np.float64)
#    vertices[:, 0] = np.cos(theta) * 0.1
#    vertices[:, 1] = np.sin(theta) * 0.1
#    vertices[:, 2] = 0.01
#
#    X, Y = np.meshgrid(x, x, indexing='ij')
#
#    points = np.zeros((3, X.size), dtype=np.float64)
#    points[0] = X.flatten()
#    points[1] = Y.flatten()
#
#    b1 = magnetic_field_current_loops(vertices, points.T, [np.arange(Ntheta)])[0]
#
#    from mayavi import mlab
#    mlab.figure()
#    q = mlab.quiver3d(*points, *b1.T)
#    q.glyph.glyph_source.glyph_position = 'center'
#
#    loops = np.array([np.arange(len(vertices)-1), np.array([3, 2, 1, 0])])
#    a1 = vector_potential(vertices, points.T, loops)[1]
#    q = mlab.quiver3d(*points, *a1.T, colormap='viridis')
#    q.glyph.glyph_source.glyph_position = 'center'
