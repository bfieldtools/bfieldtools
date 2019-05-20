import numpy as np


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy.linalg import eigh
    from mayavi import mlab
    import trimesh
    import utils
    
    from bringout_core import compute_C
    from laplacian_mesh import laplacian_matrix, mass_matrix

#    mesh = trimesh.load('./example_meshes/10x10_plane_hires.obj')
    mesh = trimesh.load('./example_meshes/disc.obj')

    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

    L = laplacian_matrix(mesh.vertices, mesh.faces)
    M = mass_matrix(mesh.vertices, mesh.faces)

#    u, v = eigh(-L.todense()[inner_verts][:,inner_verts], M.todense()[inner_verts][:,inner_verts])
    u, v = eigh(-L.todense(), M.todense())
    plt.plot(u)
#    scalars = np.zeros(L.shape[0])
#    scalars[inner_verts] = v[:, 1]
    scalars = v[:,5]
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=scalars)

    Nmodes = 30
    limit = np.max(abs(v[:,0]))

    verts= mesh.vertices
    tris = mesh.faces

    for ii in range(1,Nmodes):
        n = int(np.sqrt(Nmodes))
        i = ii % n
        j = int(ii/n)
        print(i,j)
        x = verts[:,0] + i*12
        y = verts[:,1]
        z = verts[:,2] + j*12
#        scalars[inner_verts] = v[:,ii]
#        scalars[inner] = v[:,4] +v[:,5]
        scalars = v[:,ii]
        s=mlab.triangular_mesh(x,y,z, tris, scalars=scalars, colormap = 'bwr') #M[:,70])
        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True


    
    Np = 10
    dz = 0.2
    
    x, z = np.meshgrid(np.linspace(-0.1,0.1,Np),np.linspace(-0.1,0.1,Np))
        
    fp = np.array((x.flatten(),dz*np.ones(x.flatten().shape), z.flatten())).T
    
#    mlab.figure()
#    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
#    mlab.points3d(*fp.T)
#    
    C = compute_C(mesh, fp)
    
#    vec = np.zeros(L.shape[0])
#    vec[inner_verts] = v[:, 3]
    vec = v[:,3]
    
    B = np.zeros(fp.shape)
    B[:,0] = C[:,:,0]@vec
    B[:,1] = C[:,:,1]@vec
    B[:,2] = C[:,:,2]@vec
    
        
    mlab.figure()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars = vec)
    
    mlab.quiver3d(*fp.T, *B.T)
    
    
    d = 100e-6
    sigma = 3.7e7
    T = 300
    kB = 1.38064852e-23
    mu0 = 4*np.pi*1e-7
    
    plt.figure()
    B = np.zeros(fp.shape)
    for i in range(1,5000):
#        vec = np.zeros(L.shape[0])
#        vec[inner_verts] = v[:, i]#/np.linalg.norm(v[:,i])
        vec = v[:,i]
        vec = M.sqrt()@vec
        vec *= 1/np.sqrt(u[i])
        vec *= np.sqrt(4*kB*T*sigma*d)
        B[:,0] += (C[:,:,0]@vec)**2
        B[:,1] += (C[:,:,1]@vec)**2
        B[:,2] += (C[:,:,2]@vec)**2
        
        plt.plot(i, np.sqrt(B[0,1]),'.')
        
    B = np.sqrt(B)
    
    np.mean(B, axis=0)
    
    
    Binf = mu0*np.sqrt(sigma*kB*T/(8*np.pi)*(d/(dz*(dz+d))))
    Broth = mu0*np.sqrt(kB*T*sigma*d/(8*np.pi*dz**2))