"""
Triangle potential validation
=============================
"""

import numpy as np
import matplotlib.pyplot as plt
#import sys
#path = '/m/home/home8/80/makinea1/unix/pythonstuff/bfieldtools'
#if path not in sys.path:
#    sys.path.insert(0,path)

from bfieldtools.integrals import triangle_potential_dipole_linear
from bfieldtools.integrals import omega
from bfieldtools.utils import tri_normals_and_areas

#########################################################
#%% Test potential shape slightly above the surface
points = np.array([[0,0,0],
                   [1,0,0],
                   [0,1,0]])

tris = np.array([[0,1,2]])
p_tris = points[tris]

# Evaluation points
Nx = 100
xx = np.linspace(-2, 2, Nx)
X,Y = np.meshgrid(xx, xx, indexing='ij')
Z = np.zeros_like(X) + 0.01
p_eval = np.array([X,Y,Z]).reshape(3,-1).T

# Difference vectors
RR = p_eval[:,None,None,:] - p_tris[None,:,:,:]
tn, ta = tri_normals_and_areas(points, tris)

pot = triangle_potential_dipole_linear(RR, tn, ta, False)

# Plot shapes
f, ax = plt.subplots(1, 3)
for i in range(3):
    plt.sca(ax[i])
    plt.imshow(pot[:,0,i].reshape(Nx, Nx))

#########################################################
#%% Test summation formula
pot_sum = triangle_potential_dipole_linear(RR, tn, ta, False).sum(axis=-1)
solid_angle = omega(RR)

# Plot shapes
f, ax = plt.subplots(1, 3)
plt.sca(ax[0])
plt.title('Sum of potentials')
plt.imshow(pot_sum[:,0].reshape(Nx, Nx), vmin=0, vmax=pot_sum.max())
plt.sca(ax[1])
plt.title('Solid angle')
plt.imshow(solid_angle[:,0].reshape(Nx, Nx), vmin=0, vmax=pot_sum.max())
plt.sca(ax[2])
plt.title('Difference')
plt.imshow((solid_angle[:,0]-pot_sum[:,0]).reshape(Nx, Nx),
           vmin=0, vmax=pot_sum.max())
plt.axis('image')


#########################################################
#%% Test asymptotic behavour
def dip_potential(Reval, Rdip, moment):
    R  = Reval - Rdip
    r = np.linalg.norm(R, axis=1)
    return (moment*R).sum(axis=1)/r**3

# Center of mass
Rdip = points.mean(axis=0)
# Moment
m = ta*tn
# Eval points
Neval = 100
p_eval2 = np.zeros((Neval, 3))
z = np.linspace(0,100, Neval)
p_eval2[:,2] = np.linspace(0,100, Neval)
p_eval2 += Rdip


plt.figure()

# Plot dipole field approximating uniform dipolar density
plt.plot(z, dip_potential(p_eval2, Rdip, m))
# Plot sum of the linear dipoles
RR = p_eval2[:,None,None,:] - p_tris[None,:,:,:]
pot = triangle_potential_dipole_linear(RR, tn, ta, False)
plt.plot(z,  pot.sum(axis=-1)[:,0])



