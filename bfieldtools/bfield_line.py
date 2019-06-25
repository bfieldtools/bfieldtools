import numpy as np

def mycross(r1,r2):
    """ Cross without overhead
    """
    result = np.zeros(r1.shape)
    result[0] = r1[1]*r2[2] - r1[2]*r2[1]
    result[1] = r1[2]*r2[0] - r1[0]*r2[2]
    result[2] = r1[0]*r2[1] - r1[1]*r2[0]
    return result

def bfield_line_segments(vertices, points):
    """ Compute b field of a segmented line current

        Parameters:

            vertices: (N_line, 3) vertices of the line with N_line-1 segments
            points:   (N_points, 3) Evaluation points

        Returns:

             bfield (N_points, 3) at points

        This calculation is based on integration by Griffiths
        on page 217 (3rd edition)
    """
    field = np.zeros(points.T.shape)
    for i in range(len(vertices)-1):
        r1 = vertices[i]
        r2 = vertices[i+1]
        d = ((r2-r1)/((r1-r2)**2).sum()).reshape(3,1)
        # Vectors between vertices and field points
        a1 = points.T - r1.reshape(3,1)
        a2 = points.T - r2.reshape(3,1)

        # Direction of the field
        f = mycross(a1, d)
        # Sine factor
        sinefactor = (d*a2).sum(axis=0)/np.sqrt((a2**2).sum(axis=0))
        sinefactor = sinefactor - (d*a1).sum(axis=0)/np.sqrt((a1**2).sum(axis=0))
        # Normalize direction field and divide by cylindrical distance
        s2 = (f**2).sum(axis=0)
        s2[s2 == 0] = 1e-12 # Regularize for points directly at the
                            # continuation of the line segment
        f *= (sinefactor/s2)

        field = field + f

    return field.T*1e-7

def vectorpot_line_segments(vertices, points):
    """ Compute vector potential of a segmented line current

        Returns:
            Vector potential (Npoints, Nvertices-1, 3)

    """
    segments = vertices[1:] - vertices[:-1]
    RR = points[:, None, :] - vertices[None, 1:, :]
    lengths2 = (segments**2).sum(axis=-1)
    # Normalize segments
    segments = segments/np.sqrt(lengths2)[:, None]
    # Distance vectors with directions of segments projected out
    dists = - np.cross(segments, np.cross(segments, RR))
    # Segment lengths squared over distances squared
    dists = lengths2/(dists**2).sum(axis=-1)
    return 1e-7*np.log(dists+1)[:, :, None]*segments

def flux_line_segments(vertices, vertices_loop):
    """ Compute magnetic flux created by a segmented line current
        on a another (closed loop) segmented current

        In other words, calculate mutual inductance of the currents
    """
    pass

if __name__ == "__main__":
    """ Plot field of a circular current path
    """
    x = np.linspace(-1, 1, 100)
    Ntheta = 200
    theta = np.linspace(0, 2*np.pi, Ntheta)
    vertices = np.zeros((Ntheta,3), dtype=np.float64)
    vertices[:,0] = np.cos(theta)*0.1
    vertices[:,1] = np.sin(theta)*0.1
    vertices[:,2] = 0.2

    X, Y = np.meshgrid(x, x, indexing='ij')

    points = np.zeros((3,X.size), dtype=np.float64)
    points[0] = X.flatten()
    points[1] = Y.flatten()

    b1 = bfield_line_segments(vertices, points.T)

    from mayavi import mlab
    q = mlab.quiver3d(*points, *b1.T)
    q.glyph.glyph_source.glyph_position = 'center'
