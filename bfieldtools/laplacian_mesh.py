import numpy as np
from scipy.sparse import csr_matrix, spdiags
from .utils import tri_normals_and_areas, dual_areas


def laplacian_matrix(verts, tris, tri_normals=None, tri_areas=None):
    """
    Sparse Laplace(-Beltrami) operator

    Parameters:

        verts: (n, 3) array (float)
        tris: (m, 3) array (int) - indices into the verts array


    Returns:

        Cotangent weights: w_ij = - 0.5* (cot(alpha) + cot(beta))

    """
    N = verts.shape[0]
    R = verts[tris]  # Nt x 3 (corners) x 3 (xyz)

    # Edges opposite to the vertex
    edges = np.roll(R, 1, -2) - np.roll(R, 2, -2) # Nt x 3 (edges) x 3 (x,y,z)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_normals, type(None)) or isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)
    ii = []
    jj = []
    cot = []
    # Loop over edges in triangles
    for i in range(3):
        i1 = (i+1) % 3
        i2 = (i+2) % 3
        ii.append(tris[:, i1])
        jj.append(tris[:, i2])
        cot.append(-0.5*(edges[:, i1, :]*edges[:, i2, :]).sum(axis=-1)/(2*tri_areas))

    ii = np.ravel(ii)
    jj = np.ravel(jj)
    cot = np.ravel(cot)
    # Build sparse matrix
    L = csr_matrix((cot, (ii, jj)), shape=(N, N), dtype=float)
    # Sum contribution from both triangles (alpha and beta angles)
    # neighbouring the edge
    L = L + L.T
    L = L - spdiags(L.sum(axis=0), 0, N, N)

    return L

def mass_matrix(verts, tris, tri_areas=None, da=None):
    '''
    Computes mass matrix of mesh.
    '''

    if da is None:
        if tri_areas is None:
            tri_normals, tri_areas = tri_normals_and_areas(verts, tris)

        da = dual_areas(tris, tri_areas)

    A =spdiags(da, 0, verts.shape[0], verts.shape[0]).tocsr()

    return A



if __name__ == '__main__':
    """ Example for calculating and plotting eigen modes of Laplacian
    """
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    from scipy.linalg import eigh
    from mayavi import mlab
    import os
    import trimesh

    from . import utils


    mesh = trimesh.load(os.path.join(__file__,  './example_meshes/10x10_plane_hires.obj'))

    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

    L = laplacian_matrix(mesh.vertices, mesh.faces)
    M = mass_matrix(mesh.vertices, mesh.faces)

    u, v = eigh(-L.todense()[inner_verts][:,inner_verts], M.todense()[inner_verts][:,inner_verts])

    plt.plot(u)
    scalars = np.zeros(L.shape[0])
    scalars[inner_verts] = v[:, 0]
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=scalars)

    Nmodes = 16**2
    limit = np.max(abs(v[:,0]))

    verts= mesh.vertices
    tris = mesh.faces

    for ii in range(Nmodes):
        n = int(np.sqrt(Nmodes))
        i = ii % n
        j = int(ii/n)
        print(i,j)
        x = verts[:,0] + i*12
        y = verts[:,1]
        z = verts[:,2] + j*12
        scalars[inner_verts] = v[:,ii]
#        scalars[inner] = v[:,4] +v[:,5]
        s=mlab.triangular_mesh(x,y,z, tris, scalars=scalars) #M[:,70])
        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True




#    s = mlab.triangular_mesh(*verts.T, tris, scalars=scalars)
#    s.actor.mapper.interpolate_scalars_before_mapping = True
#
#    @mlab.animate
#    def anim():
#        for i in range(100000):
#            scalars = np.zeros(L.shape[0])
#
#            prev = int(i / 10)
#            post = prev + 1
#
#            trans = (i % 10) / 10
#
#            scalars[inner_verts] = (1 - trans ) * v[:, prev]  + trans * v[:, post]
#            s.mlab_source.scalars = np.asarray(scalars, 'd')
#            yield
#
#    anim()
#    mlab.show()
