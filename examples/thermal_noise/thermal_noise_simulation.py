'''
Thermal noise computation
==========================

Three different examples:
   unit_sphere: DC Bnoise of a spherical shell at origin and comparison to analytical formula
   unit_disc: DC Bnoise of a unit disc at z-axis and comparison to analytical formula
   AC: AC Bnoise of a unit disc at one position

'''


import numpy as np
import matplotlib.pyplot as plt
import trimesh
from mayavi import mlab

from bfieldtools.mesh_properties import self_inductance_matrix, resistance_matrix
from bfieldtools.thermal_noise import compute_current_modes_ind_res, noise_covar, noise_var, visualize_current_modes
from bfieldtools.mesh_magnetics import magnetic_field_coupling

import pkg_resources


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
plt.rc('font', **font)

#Fix the simulation parameters
d = 100e-6
sigma = 3.7e7
T = 300
kB = 1.38064852e-23
mu0 = 4*np.pi*1e-7
freqs = np.array((0,))


Nchunks = 8
quad_degree = 2

##############################################################################
#Unit sphere
#------------


Np = 10
radius = np.linspace(0.1, 1, Np)
fp = np.zeros((1,3))

B = np.zeros((Np,3))
for i in range(Np):
    mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unit_sphere.stl'))
    mesh.apply_scale(radius[i])
    
    B_coupling = magnetic_field_coupling(mesh, fp, analytic = True)

    
    S = np.ones(mesh.triangles_center.shape[0])*sigma
    sheet_resistance = 1/(d*S)
    
    #Compute the resistance and inductance matrices
    R = resistance_matrix(mesh, sheet_resistance = sheet_resistance)
    M = self_inductance_matrix(mesh, Nchunks = Nchunks, quad_degree = quad_degree)
    
    vl =  compute_current_modes_ind_res(mesh,M,R, freqs, T,closed=True)
    
#    scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
#               size=(800, 800))
#    visualize_current_modes(mesh,vl[:,:,0], 8, 1)

#    vl[:,0] = np.zeros(vl[:,0].shape) # fix DC-component

    Btemp = noise_var(mesh, B_coupling, vl)
#    Btemp = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)
    B[i] = Btemp[:,:,0]

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
scene.scene.z_minus_view()
surface = scene.children[0].children[0].children[0].children[0]
surface.actor.property.representation = 'wireframe'
surface.actor.mapper.scalar_visibility = False
scene.scene.camera.position = [0.0, 0.0, -5.530686305704514]
scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [3.485379442647469, 8.118646600290083]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
scene.scene.camera.position = [0.0, 0.0, -4.570815128681416]
scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 1.0, 0.0]
scene.scene.camera.clipping_range = [2.535106977394602, 7.1443773556116374]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
mlab.savefig('/Users/joonas/Documents/Manuscripts/ThermalNoise/figures/validation/sphere.png',size=(800,800))

Ban = mu0*np.sqrt(2*sigma*d*kB*T/(3*np.pi*(radius)**2))

plt.figure(figsize = (5,5))
plt.semilogy(radius, Ban*1e15,linewidth = 2,label='Analytic')
plt.semilogy(radius, np.sqrt(B[:,2])*1e15, 'x', markersize = 10, markeredgewidth = 2, label='Numerical')
plt.grid()
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.legend(frameon = False)
plt.xlabel('Sphere radius')
plt.ylabel(r'$B_z$ noise at DC (fT/rHz)')
plt.tight_layout()


RE = np.abs((np.sqrt(B[:,2])-Ban))/np.abs(Ban)*100
plt.figure()
plt.plot(radius, RE)
plt.xlabel('Sphere radius')
plt.ylabel('Relative error (%)')

##############################################################################
#Unit disc, DC noise
#---------------------

mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unit_disc.stl'))
mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)
mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

vl = compute_current_modes(mesh)

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))

visualize_current_modes(mesh,vl, 42, 5, contours=False)

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
plt.ylabel('DC noise Bz (T/rHz)')
plt.tight_layout()

plt.figure()
plt.plot(z, np.abs((B[:,2]-Ban))/np.abs(Ban)*100)
plt.xlabel('Distance d/R')
plt.ylabel('Relative error (%)')

##############################################################################
#Closed cylinder, DC noise
#--------------------------

mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/closed_cylinder.stl'))
mesh.vertices, mesh.faces = trimesh.remesh.subdivide(mesh.vertices, mesh.faces)

    
S = np.ones(mesh.triangles_center.shape[0])*sigma
sheet_resistance = 1/(d*S)

#Compute the resistance and inductance matrices
R = resistance_matrix(mesh, sheet_resistance = sheet_resistance)
M = self_inductance_matrix(mesh, Nchunks = Nchunks, quad_degree = quad_degree)
    
vl =  compute_current_modes_ind_res(mesh,M,R, freqs, T,closed=True)

scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))

visualize_current_modes(mesh,vl[:,:,0], 8, 1)


scene = mlab.figure(None, bgcolor=(1, 1, 1), fgcolor=(0.5, 0.5, 0.5),
               size=(800, 800))
s = mlab.triangular_mesh(*mesh.vertices.T, mesh.faces)
scene.scene.z_minus_view()
surface = scene.children[0].children[0].children[0].children[0]
surface.actor.property.representation = 'wireframe'
surface.actor.mapper.scalar_visibility = False
scene.scene.isometric_view()
scene.scene.camera.position = [2.2578932293957665, 2.2578932293957665, 2.2578932293957665]
scene.scene.camera.focal_point = [0.0, 0.0, 0.0]
scene.scene.camera.view_angle = 30.0
scene.scene.camera.view_up = [0.0, 0.0, 1.0]
scene.scene.camera.clipping_range = [1.5738238620907348, 6.861972426889951]
scene.scene.camera.compute_view_plane_normal()
scene.scene.render()
mlab.savefig('/Users/joonas/Documents/Manuscripts/ThermalNoise/figures/validation/cylinder.png',size=(800,800))

Np = 30

x = np.linspace(-0.95, 0.95, Np)
fp = np.array((x,np.zeros(x.shape), np.zeros(x.shape))).T

B_coupling = magnetic_field_coupling(mesh, fp, analytic = True)
B = noise_var(mesh, B_coupling, vl)

#B = compute_dc_Bnoise(mesh,vl,fp,sigma,d,T)

a = 0.5
L = 2
rat = L/(2*a)
Gfact = 1/(8*np.pi) * ((3*rat**5+5*rat**3+2)/(rat**2*(1+rat**2)**2) + 3*np.arctan(rat))
Ban = np.sqrt(Gfact)*mu0*np.sqrt(kB*T*sigma*d)/a

plt.figure(figsize = (5,5))
plt.plot(x, Ban*np.ones(x.shape)*1e15,label='Analytic',linewidth = 2)
plt.plot(x, np.sqrt(B[:,0])*1e15,'x',label='Numerical',markersize = 10, markeredgewidth = 2,)
plt.grid()
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)
plt.legend(frameon = False)
plt.xlabel('Distance along long axis')
plt.ylabel('DC noise along axis (fT/rHz)')
plt.tight_layout()

plt.figure()
plt.semilogy(x, np.sqrt(B[:,0]),label='x')
plt.semilogy(x, np.sqrt(B[:,1]),label='y')
plt.semilogy(x, np.sqrt(B[:,2]),'--',label='z')
plt.legend()
plt.xlabel('Distance along long axis x')
plt.ylabel('DC noise (T/rHz)')



##############################################################################
#Unit disc, AC mode
#------------------

mesh = trimesh.load(pkg_resources.resource_filename('bfieldtools', 'example_meshes/unitdisc_extremelyfine.stl'))

Nfreqs = 50
freqs = np.logspace(0, 4, Nfreqs) #30 frequencies from 1 to 1000 Hz

S = np.ones(mesh.triangles_center.shape[0])*sigma
sheet_resistance = 1/(d*S)
    
#Compute the resistance and inductance matrices
R = resistance_matrix(mesh, sheet_resistance = sheet_resistance)
M = self_inductance_matrix(mesh, Nchunks = Nchunks, quad_degree = quad_degree)

vl =  compute_current_modes_ind_res(mesh,M,R, freqs, T,closed=False)

#
#fp = np.zeros((1,3))
#fp[0,2] = 0.1

Np = 10
z = np.linspace(0.05, 0.15, Np)
fp = np.array((np.zeros(z.shape), np.zeros(z.shape), z)).T

B_coupling = magnetic_field_coupling(mesh, fp, analytic = True)

Bf = np.sqrt(noise_var(mesh, B_coupling, vl))

#r = 1
#Ban = mu0*np.sqrt(sigma*d*kB*T/(8*np.pi*fp[0,2]**2))*(1/(1+fp[0,2]**2/r**2))

plt.figure()
plt.loglog(freqs,Bf[:,2,:].T*1e15,label = 'Numerical')
#plt.loglog(freqs, Ban*np.ones(freqs.shape), '--',label = 'Analytical, DC')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Field noise (fT/rHz)')
#plt.legend()
plt.grid()
#plt.grid(which='both')
plt.tight_layout()

cutf = np.zeros(Np)
for i in range(Np):
    idx = np.min(np.where(Bf[i,2,:]/Bf[i,2,0] < 1/np.sqrt(2)))
    cutf[i] = freqs[idx]

cutf_an = 1/(4*mu0*sigma*d*z)

plt.figure()
plt.plot(z, cutf_an)
plt.plot(z, cutf,'x')