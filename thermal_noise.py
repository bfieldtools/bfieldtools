import numpy as np
from scipy.linalg import eigh
from mayavi import mlab

from bringout_core import compute_C
from laplacian_mesh import laplacian_matrix, mass_matrix
from mutual_inductance_mesh import self_inductance_matrix
import utils


def compute_current_modes(mesh):
    '''
    Computes eddy-current modes for a mesh using surface laplacian.
    Uses Dirichlet boundary condition, i.e., stream function is zero at boundary:
    no current flow outside the surface.
    The modes are normalized so that the squared norm of the stream function gradient
    integrates to 1 over the surface. With this normalization, the resistances
    of the current modes are R_k = 1/(sigma*d), sigma = conductivity, d = thickness.
    
    Parameters:
        
        mesh: mesh object - the surface mesh
        
    Returns:
        
        vl: Nvertices x Nvertices array - the normalized eddy-current modes vl[:,i]
    
    '''
    
    boundary_verts, inner_verts, boundary_tris, inner_tris = utils.find_mesh_boundaries(mesh.vertices, mesh.faces, mesh.edges)
    
    L = laplacian_matrix(mesh.vertices, mesh.faces)
    M = mass_matrix(mesh.vertices, mesh.faces)
    
    u, v = eigh(-L.todense()[inner_verts][:,inner_verts], M.todense()[inner_verts][:,inner_verts])
    
    vl = np.zeros(M.shape)
    for i in range(v.shape[1]):
        vl[inner_verts,i] = v[:, i]/np.sqrt(u[i])
        
    return vl

def compute_dc_Bnoise(mesh, vl, fp, sigma, d, T):
    kB = 1.38064852e-23
    
    C = compute_C(mesh, fp)
    
    B = np.zeros(fp.shape)
    for i in range(vl.shape[1]):
        vec = vl[:,i]*np.sqrt(4*kB*T*sigma*d)
        B[:,0] += (C[:,:,0]@vec)**2
        B[:,1] += (C[:,:,1]@vec)**2
        B[:,2] += (C[:,:,2]@vec)**2
    B = np.sqrt(B)
    
    return B

def compute_ac_Bnoise(mesh, vl, fp, freqs, sigma, d, T):
    kB = 1.38064852e-23
    
    C = compute_C(mesh, fp)
    
    Mind = self_inductance_matrix(mesh.vertices,mesh.faces,mesh.face_normals)
    Mind = 0.5*(Mind+Mind.T)
    
    Mind_lap = vl.T@Mind@vl
    
    um_t, vm_t = eigh(Mind_lap)
    um_t = np.flip(um_t)
    vm_t = np.flip(vm_t,axis=1)
    
    
    um = np.zeros(um_t.shape)
    vm = np.zeros(vm_t.shape)
    
    csum = np.cumsum(um_t**2)/np.sum(um_t**2)
    Ncomps = np.max(np.where(csum < 0.99))
    
    um[0:Ncomps] = um_t[0:Ncomps]
    vm[0:Ncomps,0:Ncomps] = vm_t[0:Ncomps,0:Ncomps]
    
    B = np.zeros((freqs.shape[0], fp.shape[0],fp.shape[1]))
    
    Rk = 1/(sigma*d)
    eps = np.sqrt(4*kB*T*Rk)*np.ones((vl.shape[0]))
    
    for j in range(freqs.shape[0]):
        f = freqs[j]
        omega = 2*np.pi*f
        
        
        Rt = Rk/(Rk**2 + omega**2*um**2)
        Rt = np.diag(Rt)
        
        Lt = um/(Rk**2 + omega**2*um**2)
        Lt = np.diag(Lt)
        
        
        currents = -1*(vm@Rt@vm.T@eps + 1j*omega*vm@Lt@vm.T@eps)
        currents = np.abs(currents)
        
        for i in range(vl.shape[0]):
            vec = currents[i]*vl[:,i]
            B[j,:,0] += (C[:,:,0]@vec)**2
            B[j,:,1] += (C[:,:,1]@vec)**2
            B[j,:,2] += (C[:,:,2]@vec)**2
        print("Frequency %f computed" % (f))
        
    B = np.sqrt(B)
    
    return B

def visualize_current_modes(mesh,vl, Nmodes, scale, contours=True):
    verts= mesh.vertices
    tris = mesh.faces

    for ii in range(Nmodes):
        n = int(np.sqrt(Nmodes))
        i = ii % n
        j = int(ii/n)

        x = scale*verts[:,0] + i*12
        y = scale*verts[:,1]+ j*12
        z = scale*verts[:,2]
        
        limit = np.max(np.abs(vl[:,ii]))

        s=mlab.triangular_mesh(x,y,z, tris, scalars=vl[:,ii], colormap = 'bwr')

        s.module_manager.scalar_lut_manager.number_of_colors = 256
        s.module_manager.scalar_lut_manager.data_range = np.array([-limit,limit])
        s.actor.mapper.interpolate_scalars_before_mapping = True
        s.enable_contours = contours
        

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import trimesh
    
    example_name = 'unit_sphere'
#    example_name = 'unit_disc'
#    example_name = 'AC'
    
    d = 100e-6
    sigma = 3.7e7
    T = 300
    kB = 1.38064852e-23
    mu0 = 4*np.pi*1e-7
        
    
    if example_name == 'unit_sphere':
        Np = 10
        R = np.linspace(0.1, 1, Np)
        fp = np.zeros((1,3))
        
        B = np.zeros((Np,3))
        for i in range(Np):
            mesh = trimesh.load('./example_meshes/unit_sphere.stl')
            mesh.apply_scale(R[i])
            vl = compute_current_modes(mesh)
            
            vl[:,0] = np.zeros(vl[:,0].shape) # fix DC-component
            
            
            Btemp = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)
            B[i] = Btemp
            
            
        visualize_current_modes(mesh,vl, 40, 5)
            
        Ban = mu0*np.sqrt(2*sigma*d*kB*T/(3*np.pi*(R)**2))
        
        plt.figure()
        plt.plot(R, Ban,label='Analytic')
        plt.plot(R, B[:,2],'x',label='Numerical')
        plt.legend()
        plt.xlabel('Sphere radius')
        plt.ylabel('DC noise Bz (fT/rHz)')
        
        RE = np.abs((B[:,2]-Ban))/np.abs(Ban)*100
        plt.figure()
        plt.plot(R, np.abs((B[:,2]-Ban))/np.abs(Ban)*100)
        plt.xlabel('Sphere radius')
        plt.ylabel('Relative error (%)')
        
    elif example_name == 'unit_disc':
        mesh = trimesh.load('./example_meshes/unit_disc.stl')
        mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
        mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
        
        vl = compute_current_modes(mesh)
        
        
        visualize_current_modes(mesh,vl, 40, 5)
        
        Np = 30
        
        z = np.linspace(0.1, 1, Np)
        fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T
        
        B = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)
        
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
        
    elif example_name == 'AC':
        mesh = trimesh.load('./example_meshes/unit_disc.stl')
        mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
        mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
        
        vl = compute_current_modes(mesh)
        
        fp = np.zeros((1,3))
        fp[0,2] = 0.1
        
        Nfreqs = 30
        freqs = np.logspace(0, 3, Nfreqs)
        
        Bf = compute_ac_Bnoise(mesh,vl,fp,freqs,sigma,d,T)
        
        r = 1
        Ban = mu0*np.sqrt(sigma*d*kB*T/(8*np.pi*fp[0,2]**2))*(1/(1+fp[0,2]**2/r**2))
        
        plt.figure()
        plt.loglog(freqs,Bf[:,0,2],label = 'Numerical')
        plt.loglog(freqs, Ban*np.ones(freqs.shape), '--',label = 'Analytical, DC')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Field noise (fT/rHz)')
        plt.legend()
        plt.grid(which='both')
        plt.tight_layout()
        

    