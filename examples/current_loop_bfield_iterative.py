import numpy as np
import matplotlib.pyplot as plt
import time as t

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
    st1 = t.time()
    field = np.zeros(points.T.shape)
    for i in range(len(vertices)-1):
        r1 = vertices[i]
        r2 = vertices[i+1]
        d = ((r2-r1)/(np.square(r1-r2)).sum()).reshape(3,1)
        # Vectors between vertices and field points
        a1 = points.T - r1.reshape(3,1)
        a2 = points.T - r2.reshape(3,1)

        # Direction of the field
        f = mycross(a1, d)
        # Sine factor
        sinefactor = (d*a2).sum(axis=0)/np.sqrt(np.square(a2).sum(axis=0))
        sinefactor = sinefactor - (d*a1).sum(axis=0)/np.sqrt(np.square(a1).sum(axis=0))
        # Normalize direction field and divide by cylindrical distance
        s2 = np.square(f).sum(axis=0)
        s2[s2 == 0] = 1e-12 # Regularize for points directly at the
                            # continuation of the line segment
        f *= (sinefactor/s2)

        field = field + f
    et1 =  t.time()
    print("Execution time for line segment method is:",et1-st1)

    return field.T*1e-7



def double_factorial(n):
    if n <= 0:
        return 1
    else:
        return n * double_factorial(n-2)


def bfield_iterative(x, y, z, x_c, y_c, z_c, r, I, n):
    """ Compute b field of a current loop using an iterative method to estimate 
        elliptic integrals.
        
        Parameters:
            x, y, z: Evaluation points in 3D space. Accepts matrices and integers.
            x_c, y_c, z_c: Coordinates of the center point of the current loop.
            r: Radius of the current loop.
            I: Current of the current loop.
            n: Number of terms in the serie expansion.
            
        Returns:
            bfiels (N_points, 3) at evaluation points.
        
        This calculation is based on paper by Robert A. Schill, Jr (General 
        Relation for the Vector Magnetic Field of aCircular Current Loop: 
        A Closer Look). DOI: 10.1109/TMAG.2003.808597 
    """
    
    st2 = t.time()
    np.seterr(divide='ignore', invalid='ignore')
    u0 = 4*np.pi*1e-7
    Y = u0*I/(2*np.pi)
    
    #Change to cylideric coordinates
    rc = np.sqrt(np.power(x - x_c, 2) + np.power(y - y_c, 2))
    
    #Coefficients for estimating elliptic integrals with nth degree series 
    #expansion using Legendre polynomials
    m = 4*r*rc/(np.power((rc+r),2) + np.power((z-z_c),2))
    K = 1
    E = 1
    for i in range(1,n+1):
        K = K + np.square(double_factorial(2*i-1)/double_factorial(2*i))*np.power(m,i)
        E = E - np.square(double_factorial(2*i-1)/double_factorial(2*i))*np.power(m,i)/(2*i-1)
        
    K = K*np.pi/2
    E = E * np.pi/2

    #Calculation of radial and axial components of B-field
    Brc = Y * (z-z_c) / (rc * np.sqrt(np.power((rc+r),2) + np.power((z-z_c),2))) * (
            -K + E*(np.power(rc,2) + np.power(r,2) + np.power((z-z_c), 2)) /
            (np.power((rc-r),2) + np.power((z-z_c),2))
            )

    Bz = Y / (np.sqrt(np.power((rc+r),2) + np.power((z-z_c),2))) * (
            K - E*(np.power(rc,2) - np.power(r,2) + np.power((z-z_c), 2)) /
            (np.power((rc-r),2) + np.power((z-z_c),2))
            )
    
    #Set nan and inf values to 0
    Brc[np.isinf(Brc)] = 0 ;
    Bz[np.isnan(Bz)] = 0 ;
    Brc[np.isnan(Brc)] = 0 ;
    Bz[np.isinf(Bz)] = 0 ;
    
    #Change back to cartesian coordinates
    Bx = Brc * (x-x_c) / rc
    By = Brc * (y-y_c) / rc
    
    #Change nan values from coordinate transfer to 0
    Bx[np.isnan(Bx)] = 0;
    By[np.isnan(By)] = 0;
    
    B = np.zeros((3,X.size), dtype=np.float64)
    B[0] = Bx.flatten()
    B[1] = By.flatten()
    B[2] = Bz.flatten()
    
    et2 = t.time()
    print("Execution time for iterative method is:", et2-st2)
    
    return B.T



if __name__ == "__main__":
    """ Plot field of a circular current path
    """
    x = np.linspace(-1, 1, 100)
    Ntheta = 10000
    theta = np.linspace(0, 2*np.pi, Ntheta)
    vertices = np.zeros((Ntheta,3), dtype=np.float64)
    vertices[:,0] = np.cos(theta)*0.1
    vertices[:,1] = np.sin(theta)*0.1
    vertices[:,2] = 0.2


    X, Y = np.meshgrid(x, x, indexing='ij')
    Z = np.zeros((x.size,x.size))
    
    points = np.zeros((3,X.size), dtype=np.float64)
    points[0] = X.flatten()
    points[1] = Y.flatten()
    
    b1 = bfield_line_segments(vertices, points.T) #Calculates discretised bfield
    b2 = bfield_iterative(X, Y, Z, 0, 0, 0.2, 0.1, 1, 25) #Calculates bfield iteratively
    
    
    #Error between two calculation methods.
    berr = (b2 - b1)/b1 * 100
    BE = berr.T[2] #By changing the index, errors in different components can be obtained
    ind = np.where(np.abs(BE)>.1) #The limit for significant error is set to 0.1%
    bpoints = points.T[ind]
    
    
    
    from mayavi import mlab
    mlab.figure(1)
    q = mlab.quiver3d(*points, *b1.T)
    q.glyph.glyph_source.glyph_position = 'center'
    mlab.plot3d(*vertices.T)
    
    mlab.figure(2)
    q = mlab.quiver3d(*points, *b2.T)
    q.glyph.glyph_source.glyph_position = 'center'
    mlab.plot3d(*vertices.T)
    
    plt.figure(3)
    plt.hist(berr.T[2], bins = 50, density=True, histtype='bar')
    plt.title('Histogram of error between calculation methods.')
    plt.xlabel('%')
    
    #Plot the b-field vectors exceeding the error limit
    if len(bpoints > 0):
    
        from mayavi import mlab
        mlab.figure(3)
        q = mlab.quiver3d(*bpoints.T, *b1[ind].T)
        q.glyph.glyph_source.glyph_position = 'center'
        mlab.plot3d(*vertices.T)
    
        q = mlab.quiver3d(*bpoints.T, *b2[ind].T)
        q.glyph.glyph_source.glyph_position = 'center'

    