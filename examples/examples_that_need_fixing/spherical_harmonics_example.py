"""
Compact spherical harmonics tools and visualization example.
"""

import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt

from bfieldtools.sphtools import sphbasis, plotsph, sphfittools
from bfieldtools import sphtools


sph = sphbasis(20)

obj = plotsph.plotYlms(sph, 4)
obj = plotsph.plotYlm(sph, 3, 3)

#    obj = plotsph.plotPsilm(sph,5,2)
#    obj = plotsph.plotPhilm(sph,2,0)

offset = np.array((0, 0, 2))
mlab.figure()
obj = plotsph.plotBPhilm_volume(sph, 5, 0, 1, 10, offset)

mlab.figure()
obj = plotsph.plotBPsilm_volume(sph, 5, 0, 1, 10, offset)


Psilm1 = sphtools.Psilm(1, 0, sph.sqp[:, 1], sph.sqp[:, 2])
Psilm2 = sphtools.Psilm(7, 0, sph.sqp[:, 1], sph.sqp[:, 2])

print(sph.innerproduct(Psilm1, Psilm2))

Philm1 = sphtools.Philm(1, 0, sph.sqp[:, 1], sph.sqp[:, 2])
Philm2 = sphtools.Philm(7, 0, sph.sqp[:, 1], sph.sqp[:, 2])

print(sph.innerproduct(Philm2, Philm2))


B = np.zeros(sph.sqp.shape)
# B[:,0] = 1
B[:, 2] = sph.qp.points[:, 0] / np.max(sph.qp.points[:, 0])
B += 0.1 * np.random.randn(B.shape[0], B.shape[1])

B = sphtools.cartvec2sph(sph.sqp, B)

coeffs = sph.avsphspectra(B, 7)  # OK??

plt.figure()
plt.semilogy(coeffs ** 2)

obj = plotsph.plotYlm(sph, 5, 3)

Np = 10
lim = 3
x, y, z = np.meshgrid(
    np.linspace(-lim, lim, Np), np.linspace(-lim, lim, Np), np.linspace(-lim, lim, Np)
)


#    p = np.array((x.flatten(), y.flatten(), z.flatten())).T

p = np.array((x.flatten(), y.flatten(), np.zeros(y.flatten().shape))).T
lmax = 2
acoeffs = np.zeros(lmax * (lmax + 2))
bcoeffs = np.zeros(lmax * (lmax + 2))
acoeffs[7] = 1
#    bcoeffs[2] = 1

pot = sphtools.potential(p, acoeffs, bcoeffs, lmax)

pot = np.reshape(pot, x.shape)

mlab.mesh(x[:, :, 0], y[:, :, 0], z[:, :, 0], scalars=pot[:, :, 0], colormap="Spectral")

coords = np.zeros((p.shape[0], p.shape[1], 3))
coords[:, :, 0] = p
coords[:, :, 1] = p
coords[:, :, 2] = p

B = np.zeros((coords.shape[0], coords.shape[1]))
B[:, 2] = p[:, 0] / np.max(p[:, 0])
B[:, 1] = 0.3
B += 0.4 * np.random.randn(B.shape[0], B.shape[1])

lmax = 5
coeffs, coeffs2, mse = sphfittools.fitSpectra(sph, coords, B, lmax)

plt.figure()
plt.semilogy(coeffs ** 2, ".")
