import numpy as np


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy.linalg import eigh
    from mayavi import mlab
    import trimesh
    import utils
    
    from bringout_core import compute_C
    from laplacian_mesh import laplacian_matrix, mass_matrix
    from mutual_inductance_mesh import self_inductance_matrix
#    mesh = trimesh.load('./example_meshes/10x10_plane_hires.obj')
    mesh = trimesh.load('./example_meshes/unit_disc.stl')
#    mesh = trimesh.load('./example_meshes/unit_spiral.stl')
#    
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
#    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
#    mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
    boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)

    L = laplacian_matrix(mesh.vertices, mesh.faces)
    M = mass_matrix(mesh.vertices, mesh.faces)

    u, v = eigh(-L.todense()[inner_verts][:,inner_verts], M.todense()[inner_verts][:,inner_verts])
#    u, v = eigh(-L.todense(), M.todense())
#    plt.plot(u)
#    scalars = np.zeros(L.shape[0])
#    scalars[inner_verts] = v[:, 30]
#    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=scalars)

    Nmodes = 64
    limit = np.max(abs(v[:,0]))

    verts= mesh.vertices
    tris = mesh.faces

    for ii in range(Nmodes):
        n = int(np.sqrt(Nmodes))
        i = ii % n
        j = int(ii/n)
        print(i,j)
        x = 1*verts[:,0] + i*12
        y = 1*verts[:,1]
        z = 1*verts[:,2]+ j*12
        scalars = np.zeros(L.shape[0])
        scalars[inner_verts] = v[:,ii]
        s=mlab.triangular_mesh(x,y,z, tris, scalars=scalars, colormap = 'bwr') #M[:,70])

        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.enable_contours = True


    
    Np = 40
    dz = 0.2
#    
#    x, z = np.meshgrid(np.linspace(-0.05,0.05,Np),np.linspace(-0.05,0.05,Np))
#        
#    fp = np.array((x.flatten(),dz*np.ones(x.flatten().shape), z.flatten())).T
#    
    z = np.linspace(0.1, 10, Np)
    fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T
#    
    mlab.figure()
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
    mlab.points3d(*fp.T,scale_factor = 0.1)
    
    C = compute_C(mesh, fp)
    
#    vec = np.zeros(L.shape[0])
#    vec[inner_verts] = v[:, 3]
#    vec = v[:,3]
#    
#    B = np.zeros(fp.shape)
#    B[:,0] = C[:,:,0]@vec
#    B[:,1] = C[:,:,1]@vec
#    B[:,2] = C[:,:,2]@vec
#    
#        
#    mlab.figure()
#    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars = vec)
#    
#    mlab.quiver3d(*fp.T, *B.T)
    
    
    d = 100e-6
    sigma = 3.7e7
    T = 300
    kB = 1.38064852e-23
    mu0 = 4*np.pi*1e-7
    
    B = np.zeros(fp.shape)
    for i in range(v.shape[1]):
        vec = np.zeros(L.shape[0])
        vec[inner_verts] = v[:, i]/np.sqrt(u[i])
        vec *= np.sqrt(4*kB*T*sigma*d)
        B[:,0] += (C[:,:,0]@vec)**2
        B[:,1] += (C[:,:,1]@vec)**2
        B[:,2] += (C[:,:,2]@vec)**2
    B = np.sqrt(B)
    
#    np.mean(B, axis=0)
    
    
#    Binf = mu0*np.sqrt(sigma*kB*T/(8*np.pi)*(d/(dz*(dz+d))))
#    Broth = mu0*np.sqrt(kB*T*sigma*d/(8*np.pi*dz**2))
    r = 1
    Ban = mu0*np.sqrt(sigma*d*kB*T/(8*np.pi*z**2))*(1/(1+z**2/r**2))
    
    plt.figure()
    plt.semilogy(z, Ban,label='Analytic')
    plt.semilogy(z, B[:,2],'x',label='Numerical')
    plt.legend()
    plt.xlabel('Distance d/R')
    plt.ylabel('DC noise Bz (fT/rHz)')
    
    
    plt.figure()
    plt.plot(z, np.abs((B[:,2]-Ban))/np.abs(Ban)*100)
    plt.xlabel('Distance d/R')
    plt.ylabel('Relative error (%)')
    
    
    Mind = self_inductance_matrix(mesh.vertices,mesh.faces,mesh.face_normals)
    Mind=0.5*(Mind+Mind.T)
    vl = np.zeros(Mind.shape)
    for i in range(v.shape[1]):
        vl[inner_verts,i] = v[:, i]/np.sqrt(u[i])
    
    Mind_lap = vl@Mind@vl.T
    
    um, vm = eigh(Mind_lap)
    um = np.flip(um)
    vm = np.flip(vm,axis=1)
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=vm[:,0])
    
    
    
    
    
    M = self_inductance_matrix(mesh.vertices,mesh.faces,mesh.face_normals)
    M=0.5*(M+M.T)
    Min = M[inner_verts][:,inner_verts]
    
    L = laplacian_matrix(verts, tris)
    L = np.array(L.todense())
    
    w,v = eigh(-L[inner_verts][:,inner_verts],Min)
    
    vl = np.zeros(M.shape)
    for i in range(v.shape[1]):
        vl[inner_verts,i] = v[:, i]
#    
    mlab.triangular_mesh(*mesh.vertices.T, mesh.faces, scalars=vl[:,1])
#    
    d = 100e-6
    sigma = 3.7e7
    T = 300
    kB = 1.38064852e-23
    mu0 = 4*np.pi*1e-7
    f = 10000
    
#    eps = np.sqrt(4*kB*T/(sigma*d))*np.ones((vl.shape[0]))
    
#    A = np.linalg.pinv(1j*2*np.pi*f*Mind_lap + (1/(sigma*d))*np.eye(Mind_lap.shape[0]))
    
#    currents = A@eps
#    currents = np.abs(currents)
    
    eps = np.sqrt(4*kB*T/(sigma*d))*np.ones((vl.shape[0]))
    A = np.diag(1/np.sqrt((1/(sigma*d))**2 + (um*2*np.pi*f)**2))
    currents = vm@A@vm.T@eps
    
    B = np.zeros(fp.shape)
    for i in range(vl.shape[0]):
#        vec = np.ones(vm.shape[1]).T
#        vec *= np.sqrt(4*kB*T*sigma*d)/(np.sqrt((1 + (2*np.pi*um[i]*sigma*d*f))))
#        vec = vm@vec
#        vec = vm[:,i]
#        vec *= np.sqrt(4*kB*T*sigma*d)/(np.sqrt((1 + (2*np.pi*um[i]*sigma*d*f)**2)))
        
#        vec = vl[:,i]/w[i]
#        vec *= np.sqrt(4*kB*T/(sigma*d))
        vec = currents[i]*vl[:,i]
        B[:,0] += (C[:,:,0]@vec)**2
        B[:,1] += (C[:,:,1]@vec)**2
        B[:,2] += (C[:,:,2]@vec)**2
    B = np.sqrt(B)
    
#    np.mean(B, axis=0)
    
    
#    Binf = mu0*np.sqrt(sigma*kB*T/(8*np.pi)*(d/(dz*(dz+d))))
#    Broth = mu0*np.sqrt(kB*T*sigma*d/(8*np.pi*dz**2))
    r = 1
    Ban = mu0*np.sqrt(sigma*d*kB*T/(8*np.pi*z**2))*(1/(1+z**2/r**2))
    
    plt.figure()
    plt.plot(z, Ban,label='Analytic')
    plt.plot(z, B[:,2],'x',label='Numerical')
    plt.legend()
    plt.xlabel('Distance d/R')
    plt.ylabel('DC noise Bz (fT/rHz)')
    
    